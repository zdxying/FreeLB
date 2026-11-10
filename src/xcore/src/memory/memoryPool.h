/* This file is part of xcore
 *
 * Copyright (C) 2026 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of xcore will be updated at
 * <https://github.com/zdxying/xcore>
 *
 * xcore is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * xcore is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with xcore. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// memoryPool.h

#ifndef XCORE_MEMORY_MEMORYPOOL_H
#define XCORE_MEMORY_MEMORYPOOL_H

#include <cstddef>  // std::byte, std::size_t
#include <cstdint>  // std::uint32_t
#include <map>
#include <memory>  // std::unique_ptr
#include <set>
#include <string>
#include <vector>


namespace xcore {

// size-class scheme: power-of-two classes from 2^kClassShift to 2^(kNumClasses-1)
// bytes (8 B .. 1 MB). Requests are rounded up to the smallest class that fits
// both the size and the (rounded up) alignment.
inline constexpr std::size_t kClassShift = 3;  // smallest class = 8 bytes
inline constexpr std::size_t kNumClasses = 18;  // 8 B .. 1 MB
inline constexpr std::size_t kMaxClassSize =
  (std::size_t{1} << (kClassShift + kNumClasses - 1));  // 1 MB
// most common alignment request (operator new default); fast-path handled
// separately in allocate_raw()
inline constexpr std::size_t kDefaultAlignment = alignof(std::max_align_t);


// (cached) memoryPool for frequently used temporary objects (mainly arrays).
//
// Design overview (single-threaded, MPI/multi-process safe):
//   * small/medium sizes are served from size-class free stacks: each
//     allocation pops a fixed-size chunk (O(1)), each deallocation pushes it
//     back (O(1)). Chunks live in large contiguous super blocks (one aligned
//     operator new per block) which keeps the reused memory cache/TLB friendly.
//   * chunks are linked by an intrusive free list: the first pointer-sized
//     slot of a free chunk stores the address of the next free chunk, so no
//     separate free-list array is kept per block.
//   * chunks are assigned to the requested alignment by rounding the size
//     class up to a multiple of a power-of-two alignment; the super block base
//     is aligned to the class size, so every chunk is properly aligned.
//   * oversized requests (> kMaxClassSize) are served by a dedicated big-block
//     path; free big blocks are reused instead of being resized.
//   * fully-free super blocks are retired for reuse instead of being returned
//     to the OS (avoids alloc/free thrash); call shrink() or set an
//     auto-release threshold to give memory back on long-running workloads.
class MemoryPool {
 public:
  struct SuperBlock {
    std::byte* _base = nullptr;
    std::size_t _chunk_size = 0;
    std::uint32_t _num_chunks = 0;
    std::uint32_t _class = 0;
    std::uint32_t _free_count = 0;
    std::byte* _free_head = nullptr;  // intrusive free-list head
  };

  struct BigBlock {
    std::byte* _base = nullptr;
    std::size_t _size = 0;
    std::size_t _alignment = 0;
    bool _free = true;
  };

  struct Region {
    std::byte* _base = nullptr;
    std::size_t _size = 0;
    SuperBlock* _sb = nullptr;  // exactly one of _sb / _big is set
    BigBlock* _big = nullptr;
  };

  MemoryPool(std::size_t est_size = 128) {
    _active.resize(kNumClasses);
    _retired.resize(kNumClasses);
    _regions.reserve(est_size);
    _blocks.reserve(est_size);
  }

  MemoryPool(const MemoryPool&) = delete;
  MemoryPool& operator=(const MemoryPool&) = delete;
  ~MemoryPool();

  void* allocate_raw(std::size_t size,
                     std::size_t alignment = alignof(std::max_align_t));

  void deallocate_raw(void* ptr);

  void print_status(int level = 0) const;

  // return all fully-free super blocks and free big blocks to the OS
  void shrink();

  // if the pool holds more than `bytes` while freeing, release memory
  // immediately instead of retiring it (0 disables auto-release, default)
  void set_auto_release_threshold(std::size_t bytes) {
    _auto_release_threshold = bytes;
  }

  std::size_t getTotalSize() const { return _total_size; }
  std::size_t getUsedSize() const { return _used_size; }

  static MemoryPool& getInstance() {
    static MemoryPool instance;
    return instance;
  }

 private:
  std::byte* allocate_class(std::size_t eff);
  std::byte* allocate_big(std::size_t size, std::size_t alignment);
  SuperBlock* add_superblock(std::size_t idx);
  void release_superblock(SuperBlock* sb);
  std::size_t find_region(std::byte* ptr) const;
  void register_region(std::byte* base, std::size_t size, void* owner,
                       bool is_big);
  void unregister_region(std::byte* base);

  // super blocks (chunk pools) and big blocks; stable object addresses, so
  // _active/_retired/_regions can hold raw pointers
  std::vector<std::unique_ptr<SuperBlock>> _blocks;
  std::vector<std::unique_ptr<BigBlock>> _big_blocks;
  // stack of pointers to currently-free big blocks (for O(#free) reuse scan)
  std::vector<BigBlock*> _big_free;
  // per-class stacks of super blocks that still have free chunks (_active) or
  // are completely free and waiting for reuse (_retired)
  std::vector<std::vector<SuperBlock*>> _active;
  std::vector<std::vector<SuperBlock*>> _retired;
  // sorted by _base; maps a pointer back to its owner block
  std::vector<Region> _regions;

  std::size_t _total_size = 0;
  std::size_t _used_size = 0;
  std::size_t _high_water = 0;
  std::size_t _auto_release_threshold = 0;
  // index into _regions of the last block touched by deallocate_raw(); let
  // the common roll/reuse pattern skip the region binary search
  std::uint32_t _hint = 0;
};

// memoryPoolTracker for tracking MemoryPool's allocations and deallocations
class MemoryPoolTracker {
 public:
  static MemoryPoolTracker& getInstance() {
    static MemoryPoolTracker instance;
    return instance;
  }

  void on_allocate(std::size_t size);
  void on_deallocate(std::size_t size);

  class Scope {
   public:
    Scope(const std::string& rec_name = "default");
    ~Scope();

    void start();
    void stop();
    void print(int level = 0) const;

   private:
    std::string _rec_name;
    MemoryPoolTracker& _tracker;
  };

 private:
  MemoryPoolTracker() = default;
  MemoryPoolTracker(const MemoryPoolTracker&) = delete;
  MemoryPoolTracker& operator=(const MemoryPoolTracker&) = delete;

  std::map<std::string, std::vector<std::size_t>> _alloc_rec;
  std::map<std::string, std::vector<std::size_t>> _dealloc_rec;
  std::set<std::string> _active_records;
};

}  // namespace xcore
#endif  // XCORE_MEMORY_MEMORYPOOL_H
