# xcore MemoryPool Design

`xcore::MemoryPool` is a high-performance, single-threaded pool allocator aimed at
frequently created temporary objects (especially arrays). It trades a bounded
amount of retained memory for O(1) allocation/deallocation and cache-friendly
reuse. This document explains its core design.

## Size Classes

Requests up to 1 MB are rounded up to one of 18 power-of-two size classes
(`kNumClasses = 18`), starting at 8 bytes (`kClassShift = 3`). A request of `n`
bytes maps to the smallest class `2^k` with `2^k >= n` via `to_class()`
(essentially `ceil_log2(n)`). Rounding up means a class stores chunks of a fixed
size, so no per-chunk metadata or free-list header is needed.

## Super Blocks

Each class's memory lives in **super blocks** — large contiguous ranges acquired
with a single aligned `operator new`. A super block is divided into
`num_chunks` fixed-size chunks, where `num_chunks` is chosen to target ~64 KB of
total memory per block (`kSuperBlockTarget`), clamped to `[1, 512]`.

Because every chunk is identical in size and the block base is aligned to the
chunk size, free chunks are linked by an **intrusive free list**: the first
pointer-sized slot of a free chunk stores the address of the next free chunk, and
the block keeps only a head pointer (`_free_head`) plus a free count. No separate
free-list array is allocated per block — the list lives inside the chunks
themselves, which keeps overhead near zero and the reused memory cache friendly.
The blocks are cache/TLB friendly since reuse stays within few large regions.

## Allocation Flow

`allocate_raw(size, alignment)`:

1. Reject `size == 0`; normalize a zero alignment to `kDefaultAlignment`.
2. If `alignment == kDefaultAlignment` (common case), set `eff = max(size, kDefaultAlignment)`. Otherwise round the alignment up to a power of two `align_pow` and set `eff = max(size, align_pow)`.
3. If `eff <= kMaxClassSize`, go through the size-class path (`allocate_class`), else the big-block path (`allocate_big`).

`allocate_class(eff)`:

1. Map `eff` to its class index `idx`.
2. If no active block exists for `idx`, pop a retired super block of that class, or if none, create one with `add_superblock`. Push it onto `_active[idx]`.
3. Pop the head chunk from the active block's intrusive free list: `p = free_head; free_head = *(void**)p`. If the block is now exhausted (`--free_count == 0`), remove it from `_active[idx]`.
4. Account `_used_size += chunk_size`.

`add_superblock(idx)` allocates a new `SuperBlock` object, then `chunk * num`
bytes with a single aligned `operator new` (`align_val_t{chunk}`). It threads the
intrusive free list through all chunks (chunk *i* points to chunk *i+1*) and
registers the region. `_blocks` holds `std::unique_ptr<SuperBlock>` so object
addresses stay stable and released blocks free themselves.

`allocate_big(size, alignment)` scans `_big_free` for a previously freed block
with `_size >= size` and `_alignment >= alignment` and reuses it; otherwise it
allocates a fresh aligned block and registers it.

## Deallocation Flow

`deallocate_raw(ptr)`:

1. Ignore `nullptr`.
2. Locate the owning region: first test the `_hint` region (cheap pointer-range
   check), and only fall back to the sorted-array binary search `find_region` on
   a miss. An unmatched pointer is reported and dropped.
3. **Big block:** subtract its `_size` from `_used_size`. If the auto-release
   threshold is exceeded, free it to the OS immediately; otherwise mark it free
   and push its pointer onto `_big_free`.
4. **Super block:** verify `(ptr - base) % chunk_size == 0` (a chunk-aligned
   address), push the chunk onto the intrusive free list
   (`*(void**)p = free_head; free_head = p`), and subtract `chunk_size` from
   `_used_size`.
5. If the block is now completely free (`++free_count == num_chunks`), remove it
   from `_active[idx]`; then either release it to the OS (threshold exceeded) or
   move it to `_retired[idx]` for reuse.

Both core operations are O(1) amortized: a pop/push on the intrusive free list
plus a cached (often skipped) region lookup.

## Alignment

Alignment is folded into the class selection rather than handled with headers.
An over-aligned request rounds `alignment` up to a power of two `align_pow`, then
uses the class `max(size, align_pow)`. Since the super block base is aligned to
its chunk size, every chunk in the block is aligned to at least `align_pow`. The
common default-alignment (`alignof(std::max_align_t)`) path is handled separately
as a fast path.

## Big Blocks

Requests larger than 1 MB (`kMaxClassSize`) take a dedicated big-block path: a
single aligned `operator new` per request. Freed big blocks are kept on a free
stack (`_big_free`) and reused for later requests of equal-or-smaller size and
alignment, avoiding repeated `mmap`/`brk` churn.

## Memory Reclamation

Retired super blocks and free big blocks can be returned to the OS with
`shrink()`. For long-running processes, `set_auto_release_threshold(bytes)`
makes the pool release a fully free block immediately once `_total_size` exceeds
the threshold, bounding peak retention.

## Region Lookup

`deallocate_raw()` must map a raw pointer back to its owning block. All allocated
regions (super blocks and big blocks) are kept sorted by base address in
`_regions`, enabling a binary search (`find_region`). Because a typical workload
frees in a roll/reuse order, a `_hint` index caches the last-touched region and
skips the binary search on the common path.

A freed pointer is validated: big-block ownership and chunk-boundary alignment
(`off % chunk_size == 0`) are checked, and out-of-pool pointers are reported via
`ERROR_MESSAGE` instead of corrupting the heap.

Blocks are owned as stable heap objects (`std::vector<std::unique_ptr<SuperBlock>>`
and the big-block counterpart), so `_active`/`_retired`/`_regions` store raw
pointers directly — no index arithmetic or released-slot freelists are needed.

## Usage

- Singleton access via `MemoryPool::getInstance()`, or instantiate a private pool.
- `allocate_raw(size, alignment)` / `deallocate_raw(ptr)` for raw allocation.
- `PoolAllocator<T>` adapts the pool to the standard allocator interface, so it
  can be plugged into `std::vector`, `std::list`, etc.
- `MemoryPoolTracker` (enabled by `MEMPOOL_TRACK`) records per-scope allocation
  totals for profiling.
- `print_status()` reports used/total/high-water bytes for tuning.

## Design Trade-offs

- **Single-threaded.** No locking or per-thread caching; safe for MPI/multi-process
  use where each process has its own pool.
- **Retains memory.** Fully free blocks are held for reuse until `shrink()` or an
  auto-release threshold; ideal for repeated alloc/free workloads, not for the
  coldest memory profile.
- **Rounding waste.** Small allocations pay the size-class rounding (e.g. 10 bytes
  costs a 16-byte chunk), typical of pool allocators.
