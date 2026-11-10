# xcore MemoryPool Design

Status: current as of 2026-10-31. Exit: obsolete when `MemoryPool` is replaced or its allocation/reclamation policy changes.

`xcore::MemoryPool` is a single-threaded pool allocator for frequently created
temporary objects; it trades bounded retained memory for fast, cache-friendly reuse.

## Size Classes

Requests up to 1 MB are rounded up to one of 18 power-of-two classes
(`kNumClasses = 18`) from 8 bytes (`kClassShift = 3`): `to_class()` returns the
smallest `2^k >= n`. All chunks in a class share one size, so the free-list link
lives inside the free chunk; no per-chunk metadata is stored outside the chunk
(`memoryPool.cpp`, `add_superblock`).

## Super Blocks

Each class owns **super blocks**: one aligned `operator new` per contiguous
range, split into equal chunks whose count targets ~64 KB (`kSuperBlockTarget`,
clamped to `[1, 512]`). Free chunks form an **intrusive free list**: `_free_head`
plus a count per block, each free chunk's first pointer-sized slot holding the
next. No separate free-list array, so overhead is low and reuse is cache/TLB friendly.

## Allocation Flow

`allocate_raw()` (`memoryPool.cpp`) computes `eff = max(size, alignment rounded
up to a power of two)`: `eff <= kMaxClassSize` takes `allocate_class`, else
`allocate_big`. Class allocation pops the active free-list head, reuses a retired
block before creating one, and adds the chunk size to `_used_size`; big
allocation reuses a fitting `_big_free` entry (linear scan) or allocates fresh.

## Deallocation Flow

`deallocate_raw()` resolves the region through `_hint`, then binary-search
`find_region`; an out-of-region pointer or a super-block pointer not on a chunk
boundary is fatal: `ERROR_MESSAGE` is `[[noreturn]]` and calls `std::abort()`
(`util/assert.h`), terminating the process. A big block is released to the OS
above the auto-release threshold, else pushed on `_big_free`; a super-block chunk
is pushed on the free list, and a fully free block is released or retired. The
small-object fast path is O(1) amortized: free-list push/pop plus a usually-skipped
`_hint` check (a miss costs O(log #regions)); it is not O(1) in general because
`allocate_big` scans `_big_free` (O(#free)) and retiring a super block scans the
class's active list. No memory-pool benchmark exists in this tree, so these are
code-level claims (`memoryPool.cpp`: `allocate_big`, `deallocate_raw`).

## Alignment

Alignment is folded into class selection, not handled with headers: an
over-aligned request rounds `alignment` up to `align_pow` and uses class
`max(size, align_pow)`; the chunk-aligned base makes every chunk meet it, and
default alignment (`alignof(std::max_align_t)`) is a separate fast path.

## Big Blocks

Requests with `max(size, rounded-up alignment) > 1 MB` (`kMaxClassSize`) take the
big-block path: one aligned `operator new` per request, including a small
over-aligned request whose alignment, not size, exceeds 1 MB
(`memoryPool.cpp:94-98`). Freed blocks stay on `_big_free` for reuse at
equal-or-smaller size and alignment, avoiding repeated `mmap`/`brk` churn.

## Memory Reclamation

`shrink()` returns retired super blocks and free big blocks to the OS.
`set_auto_release_threshold(bytes)` releases a fully free block immediately once
`_total_size` exceeds the threshold, bounding peak retention.

## Region Lookup

`_regions` is sorted by base, so `find_region` binary-searches; `_hint` caches the
last-touched region. Stable `std::unique_ptr` block objects let the lists hold raw pointers.

## Usage

- Singleton access via `MemoryPool::getInstance()`, or instantiate a private pool.
- `allocate_raw(size, alignment)` / `deallocate_raw(ptr)` for raw allocation.
- `PoolAllocator<T>` adapts the pool to the standard allocator interface, so it
  can be plugged into `std::vector`, `std::list`, etc.
- `MemoryPoolTracker` records per-scope totals only when `MEMPOOL_TRACK` is
  defined; xcore defines it only for `DEBUG=1` (`src/xcore/Makefile:32`), and
  FreeLB links the release library (`src/xcore/build/install`, `make.mk`), so
  tracking is off by default.
- `print_status()` reports used/total/high-water bytes for tuning.

## Design Trade-offs

- **Single-threaded.** No locking or per-thread caching; safe for MPI/multi-process
  use with one pool per process.
- **Retains memory.** Fully free blocks are held until `shrink()` or an
  auto-release threshold; good for repeated alloc/free, not cold workloads.
- **Rounding waste.** Small allocations pay class rounding (10 bytes cost a 16-byte chunk).

## Rejected Alternatives

- **Return fully free super blocks to the OS immediately.** Rejected: repeated
  alloc/free would thrash `mmap`/`munmap`. Blocks are retired for reuse instead;
  `shrink()` or `set_auto_release_threshold()` makes release explicit. Evidence:
  `_retired` in `memoryPool.h` and the `release_superblock()` branch in
  `deallocate_raw()`.
