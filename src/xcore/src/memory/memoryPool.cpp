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

// memoryPool.cpp
#include "memory/memoryPool.h"

#include <algorithm>  // std::clamp, std::max
#include <iostream>
#include <new>  // std::align_val_t

#include "util/assert.h"


namespace xcore {

static constexpr std::size_t kb = 1024;
static constexpr std::size_t mb = kb * kb;
static constexpr std::size_t gb = mb * kb;

// one aligned operator new per super block, targeting ~64 KB of chunks
static constexpr std::size_t kSuperBlockTarget = std::size_t{64} << 10;
static constexpr std::size_t kMinChunksPerBlock = 1;
static constexpr std::size_t kMaxChunksPerBlock = 512;

namespace {

#if defined(__GNUC__) || defined(__clang__)
std::size_t ceil_log2(std::size_t x) {
  return x <= 1 ? 0 : std::size_t{64} - static_cast<std::size_t>(__builtin_clzll(x - 1));
}
#else
std::size_t ceil_log2(std::size_t x) {
  std::size_t r = 0;
  while ((std::size_t{1} << r) < x) ++r;
  return r;
}
#endif

std::size_t class_size(std::size_t idx) { return std::size_t{1} << (idx + kClassShift); }

// smallest class that holds `size` bytes (8 B .. 1 MB)
std::size_t to_class(std::size_t size) {
  if (size <= class_size(0)) return 0;
  const std::size_t idx = ceil_log2(size) - kClassShift;
  return idx < kNumClasses ? idx : kNumClasses - 1;
}

}  // namespace


MemoryPool::~MemoryPool() {
  for (auto& b : _blocks) {
    if (b->_base != nullptr) {
      ::operator delete(b->_base, std::align_val_t{b->_chunk_size});
    }
  }
  for (auto& b : _big_blocks) {
    if (b->_base != nullptr) {
      ::operator delete(b->_base, std::align_val_t{b->_alignment});
    }
  }
}

void* MemoryPool::allocate_raw(std::size_t size, std::size_t alignment) {
  if (size == 0) return nullptr;
  const std::size_t a = (alignment == 0) ? kDefaultAlignment : alignment;
#ifdef MEMPOOL_TRACK
  MemoryPoolTracker::getInstance().on_allocate(size);
#endif

  if (a == kDefaultAlignment) {
    // common fast path: default alignment, plain size
    const std::size_t eff = std::max(size, kDefaultAlignment);
    if (eff <= kMaxClassSize) return allocate_class(eff);
    return allocate_big(size, kDefaultAlignment);
  }
  // over-aligned request: fold alignment into the class size
  const std::size_t align_pow = std::size_t{1} << ceil_log2(a);
  const std::size_t eff = std::max(size, align_pow);
  if (eff <= kMaxClassSize) return allocate_class(eff);
  return allocate_big(size, a);
}

std::byte* MemoryPool::allocate_class(std::size_t eff) {
  const std::size_t idx = to_class(eff);
  // grab (or create) a super block of this class with free chunks
  if (_active[idx].empty()) {
    SuperBlock* sb;
    if (!_retired[idx].empty()) {
      sb = _retired[idx].back();
      _retired[idx].pop_back();
    } else {
      sb = add_superblock(idx);
    }
    _active[idx].push_back(sb);
  }
  SuperBlock* sb = _active[idx].back();
  std::byte* ptr = sb->_free_head;
  sb->_free_head = *reinterpret_cast<std::byte**>(ptr);
  if (--sb->_free_count == 0) _active[idx].pop_back();
  _used_size += sb->_chunk_size;
  return ptr;
}

std::byte* MemoryPool::allocate_big(std::size_t size, std::size_t alignment) {
  // oversized request: big-block path (free blocks are reused, never resized)
  for (std::size_t fi = 0; fi < _big_free.size(); ++fi) {
    BigBlock* b = _big_free[fi];
    if (b->_base != nullptr && b->_size >= size && b->_alignment >= alignment) {
      _big_free[fi] = _big_free.back();
      _big_free.pop_back();
      b->_free = false;
      _used_size += b->_size;
      return b->_base;
    }
  }
  auto bb = std::make_unique<BigBlock>();
  BigBlock* b = bb.get();
  b->_base = static_cast<std::byte*>(
    ::operator new(size, std::align_val_t{alignment}));
  b->_size = size;
  b->_alignment = alignment;
  b->_free = false;
  register_region(b->_base, b->_size, b, true);
  _big_blocks.push_back(std::move(bb));
  _total_size += size;
  _high_water = std::max(_high_water, _total_size);
  _used_size += size;
  return b->_base;
}

void MemoryPool::deallocate_raw(void* ptr) {
  if (ptr == nullptr) return;
  std::byte* p = static_cast<std::byte*>(ptr);

  // fast path: the pointer is usually in the block freed just before it
  std::size_t r = _regions.size();
  if (_hint < _regions.size()) {
    const Region& h = _regions[_hint];
    if (p >= h._base && p < h._base + h._size) r = _hint;
  }
  if (r == _regions.size()) {
    r = find_region(p);
    if (r == _regions.size()) {
      ERROR_MESSAGE("[MemoryPool] Error: Deallocating a non-existing pointer");
      return;
    }
  }
  _hint = static_cast<std::uint32_t>(r);
  const Region& reg = _regions[r];

  if (reg._big != nullptr) {
    BigBlock* b = reg._big;
    _used_size -= b->_size;
#ifdef MEMPOOL_TRACK
    MemoryPoolTracker::getInstance().on_deallocate(b->_size);
#endif
    if (_auto_release_threshold != 0 && _total_size > _auto_release_threshold &&
        b->_free == false) {
      unregister_region(b->_base);
      ::operator delete(b->_base, std::align_val_t{b->_alignment});
      _total_size -= b->_size;
      b->_base = nullptr;
      b->_free = true;
      return;
    }
    b->_free = true;
    _big_free.push_back(b);
    return;
  }

  SuperBlock* sb = reg._sb;
  const std::size_t off = static_cast<std::size_t>(p - sb->_base);
  if (sb->_base == nullptr || off % sb->_chunk_size != 0) {
    ERROR_MESSAGE("[MemoryPool] Error: Deallocating a non-existing pointer");
    return;
  }
  *reinterpret_cast<std::byte**>(p) = sb->_free_head;
  sb->_free_head = p;
  _used_size -= sb->_chunk_size;
#ifdef MEMPOOL_TRACK
  MemoryPoolTracker::getInstance().on_deallocate(sb->_chunk_size);
#endif

  // super block became completely free: retire it for reuse
  if (++sb->_free_count == sb->_num_chunks) {
    auto& act = _active[sb->_class];
    auto it = std::find(act.begin(), act.end(), sb);
    if (it != act.end()) act.erase(it);

    if (_auto_release_threshold != 0 && _total_size > _auto_release_threshold) {
      release_superblock(sb);
    } else {
      _retired[sb->_class].push_back(sb);
    }
  }
}

void MemoryPool::print_status(int level) const {
  std::cout << "[MemoryPool] Status:\n";

  std::size_t n_blocks = 0;
  std::size_t n_big = 0;
  for (const auto& b : _blocks) {
    if (b->_base != nullptr) ++n_blocks;
  }
  for (const auto& b : _big_blocks) {
    if (b->_base != nullptr) ++n_big;
  }
  std::cout << "  SuperBlocks: " << n_blocks << " | BigBlocks: " << n_big << "\n";

  if (level > 0) {
    for (std::size_t idx = 0; idx < kNumClasses; ++idx) {
      const std::size_t chunk = class_size(idx);
      std::size_t n_active = 0;
      std::size_t n_free_chunks = 0;
      for (SuperBlock* sb : _active[idx]) {
        n_active += 1;
        n_free_chunks += sb->_free_count;
      }
      for (SuperBlock* sb : _retired[idx]) {
        n_free_chunks += sb->_num_chunks;
      }
      if (n_active || n_free_chunks) {
        std::cout << "    class " << chunk << " B: " << n_active << " active, "
                  << n_free_chunks << " free chunks\n";
      }
    }
  }

  std::cout << "  Used " << _used_size << " bytes / Total " << _total_size << " bytes";
  const std::size_t total_size = _total_size;
  if (total_size >= kb) {
    std::cout << " | " << double(total_size) / kb << " KB";
    if (total_size >= mb) {
      std::cout << " | " << double(total_size) / mb << " MB";
      if (total_size >= gb) {
        std::cout << " | " << double(total_size) / gb << " GB";
      }
    }
  }
  std::cout << "\n";
  std::cout << "  High water: " << _high_water << " bytes\n";
  std::cout << std::endl;
}

void MemoryPool::shrink() {
  for (std::size_t idx = 0; idx < kNumClasses; ++idx) {
    for (SuperBlock* sb : _retired[idx]) release_superblock(sb);
    _retired[idx].clear();
  }
  const std::vector<BigBlock*> to_release = _big_free;
  for (BigBlock* b : to_release) {
    if (b->_base != nullptr) {
      unregister_region(b->_base);
      ::operator delete(b->_base, std::align_val_t{b->_alignment});
      _total_size -= b->_size;
      b->_base = nullptr;
      b->_free = true;
    }
  }
  _big_free.clear();
}

MemoryPool::SuperBlock* MemoryPool::add_superblock(std::size_t idx) {
  auto block = std::make_unique<SuperBlock>();
  SuperBlock* sb = block.get();
  const std::size_t chunk = class_size(idx);
  const std::uint32_t num = static_cast<std::uint32_t>(std::clamp(
    kSuperBlockTarget / chunk, kMinChunksPerBlock, kMaxChunksPerBlock));
  sb->_chunk_size = chunk;
  sb->_num_chunks = num;
  sb->_class = static_cast<std::uint32_t>(idx);
  sb->_free_count = num;
  sb->_base = static_cast<std::byte*>(
    ::operator new(chunk * num, std::align_val_t{chunk}));
  // build the intrusive free list: chunk i points to chunk i+1
  sb->_free_head = sb->_base;
  for (std::uint32_t i = 0; i < num; ++i) {
    std::byte* slot = sb->_base + i * chunk;
    *reinterpret_cast<std::byte**>(slot) =
      (i + 1 < num) ? slot + chunk : nullptr;
  }
  register_region(sb->_base, chunk * num, sb, false);
  _blocks.push_back(std::move(block));
  _total_size += chunk * num;
  _high_water = std::max(_high_water, _total_size);
  return sb;
}

void MemoryPool::release_superblock(SuperBlock* sb) {
  if (sb == nullptr || sb->_base == nullptr) return;
  unregister_region(sb->_base);
  ::operator delete(sb->_base, std::align_val_t{sb->_chunk_size});
  _total_size -= sb->_chunk_size * sb->_num_chunks;
  sb->_base = nullptr;
  sb->_free_count = 0;
  sb->_free_head = nullptr;
  auto it = std::find_if(_blocks.begin(), _blocks.end(),
                         [sb](const std::unique_ptr<SuperBlock>& p) {
                           return p.get() == sb;
                         });
  if (it != _blocks.end()) _blocks.erase(it);
}

std::size_t MemoryPool::find_region(std::byte* ptr) const {
  // binary search for the region whose [base, base+size) range holds `ptr`
  std::size_t lo = 0, hi = _regions.size();
  while (lo < hi) {
    const std::size_t mid = (lo + hi) / 2;
    if (ptr < _regions[mid]._base) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  if (lo == 0) return _regions.size();
  const Region& r = _regions[lo - 1];
  if (ptr < r._base || ptr >= r._base + r._size) return _regions.size();
  return lo - 1;
}

void MemoryPool::register_region(std::byte* base, std::size_t size, void* owner,
                                 bool is_big) {
  auto it = std::lower_bound(
    _regions.begin(), _regions.end(), base,
    [](const Region& r, const std::byte* b) { return r._base < b; });
  Region reg;
  reg._base = base;
  reg._size = size;
  if (is_big) {
    reg._big = static_cast<BigBlock*>(owner);
  } else {
    reg._sb = static_cast<SuperBlock*>(owner);
  }
  _regions.insert(it, reg);
}

void MemoryPool::unregister_region(std::byte* base) {
  auto it = std::lower_bound(
    _regions.begin(), _regions.end(), base,
    [](const Region& r, const std::byte* b) { return r._base < b; });
  if (it != _regions.end() && it->_base == base) {
    const std::size_t idx = static_cast<std::size_t>(it - _regions.begin());
    _regions.erase(it);
    // indices at/after `idx` shifted; invalidate the hint to stay safe
    if (_hint >= idx) _hint = static_cast<std::uint32_t>(_regions.size());
  }
}

// ---------------------------------------------------------------------------
// MemoryPoolTracker
// ---------------------------------------------------------------------------

void MemoryPoolTracker::on_allocate(std::size_t size) {
  for (const auto& rec_name : _active_records) {
    _alloc_rec[rec_name].push_back(size);
  }
}

void MemoryPoolTracker::on_deallocate(std::size_t size) {
  for (const auto& rec_name : _active_records) {
    _dealloc_rec[rec_name].push_back(size);
  }
}

MemoryPoolTracker::Scope::Scope(const std::string& rec_name)
    : _rec_name(rec_name), _tracker(MemoryPoolTracker::getInstance()) {
  _tracker._active_records.insert(rec_name);
  _tracker._alloc_rec[rec_name].clear();
  _tracker._dealloc_rec[rec_name].clear();
}

MemoryPoolTracker::Scope::~Scope() {
  stop();
  _tracker._alloc_rec.erase(_rec_name);
  _tracker._dealloc_rec.erase(_rec_name);
}

void MemoryPoolTracker::Scope::start() {
  _tracker._active_records.insert(_rec_name);
  _tracker._alloc_rec[_rec_name].clear();
  _tracker._dealloc_rec[_rec_name].clear();
}

void MemoryPoolTracker::Scope::stop() { _tracker._active_records.erase(_rec_name); }

void MemoryPoolTracker::Scope::print(int level) const {
  if (_tracker._alloc_rec.find(_rec_name) == _tracker._alloc_rec.end()) {
    std::cout << "[MemoryPoolTracker]: name " << _rec_name << " not found" << "\n";
    return;
  }
  std::cout << "[MemoryPoolTracker] " << _rec_name << ":\n";
  std::size_t total_alloc{}, total_dealloc{};
  const auto& alloc_rec = _tracker._alloc_rec[_rec_name];
  for (auto size : alloc_rec) {
    total_alloc += size;
  }
  const auto& dealloc_rec = _tracker._dealloc_rec[_rec_name];
  for (auto size : dealloc_rec) {
    total_dealloc += size;
  }

  std::cout << "  Allocated " << total_alloc << " bytes";
  if (total_alloc >= kb) {
    std::cout << " | " << double(total_alloc) / kb << " KB";
    if (total_alloc >= mb) {
      std::cout << " | " << double(total_alloc) / mb << " MB";
      if (total_alloc >= gb) {
        std::cout << " | " << double(total_alloc) / gb << " GB";
      }
    }
  }
  std::cout << "\n";
  std::cout << "  Allocation count: " << alloc_rec.size() << "\n";
  if (level > 0) {
    std::cout << "  Allocation records(in bytes):\n";
    for (std::size_t i = 0; i < alloc_rec.size(); ++i) {
      std::cout << "    " << i << ": " << alloc_rec[i] << "\n";
    }
  }

  std::cout << "  Deallocated " << total_dealloc << " bytes";
  if (total_dealloc >= kb) {
    std::cout << " | " << double(total_dealloc) / kb << " KB";
    if (total_dealloc >= mb) {
      std::cout << " | " << double(total_dealloc) / mb << " MB";
      if (total_dealloc >= gb) {
        std::cout << " | " << double(total_dealloc) / gb << " GB";
      }
    }
  }
  std::cout << "\n";
  std::cout << "  Deallocation count: " << dealloc_rec.size() << "\n";
  if (level > 0) {
    std::cout << "  Deallocation records(in bytes):\n";
    for (std::size_t i = 0; i < dealloc_rec.size(); ++i) {
      std::cout << "    " << i << ": " << dealloc_rec[i] << "\n";
    }
  }

  std::cout << std::endl;
}

}  // namespace xcore