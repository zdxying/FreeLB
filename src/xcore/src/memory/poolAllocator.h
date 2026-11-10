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

// poolAllocator.h

#ifndef XCORE_MEMORY_POOLALLOCATOR_H
#define XCORE_MEMORY_POOLALLOCATOR_H

#include <new>  // placement new

#include "xcore/src/memory/memoryPool.h"


namespace xcore {

template <typename T>
class PoolAllocator {
 public:
  // compatible with c++'s allocator
  using value_type = T;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  template <typename U>
  struct rebind {
    using other = PoolAllocator<U>;
  };

  PoolAllocator() = default;

  template <typename U>
  PoolAllocator(const PoolAllocator<U>&) noexcept {}

  constexpr size_type max_size() noexcept { return std::size_t(-1) / sizeof(T); }

  friend bool operator==(const PoolAllocator&, const PoolAllocator&) noexcept {
    return true;
  }

  friend bool operator!=(const PoolAllocator&, const PoolAllocator&) noexcept {
    return false;
  }

  using pool = MemoryPool;
  using type = PoolAllocator<T>;

  T* allocate(std::size_t n) {
    return static_cast<T*>(pool::getInstance().allocate_raw(n * sizeof(T)));
  }

  void deallocate(T* ptr, std::size_t) noexcept {
    pool::getInstance().deallocate_raw(ptr);
  }

  template <typename U, typename... Args>
  void construct(U* p, Args&&... args) {
    ::new ((void*)p) U(std::forward<Args>(args)...);
  }

  template <typename U>
  void destroy(U* p) noexcept { p->~U(); }


 private:
};

}  // namespace xcore
#endif  // XCORE_MEMORY_POOLALLOCATOR_H