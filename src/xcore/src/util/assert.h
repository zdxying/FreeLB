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

// assert.h

#ifndef XCORE_UTIL_ASSERT_H
#define XCORE_UTIL_ASSERT_H

#include <stdlib.h>

#include <iostream>
#include <sstream>

namespace xcore {

namespace detail {

[[noreturn]] inline void assert_fail_impl(
  const char* expr, const char* file, int line, const char* func) {
  std::cerr << "Assertion failed: " << expr << '\n'
            << "File: " << file << '\n'
            << "Line: " << line << '\n'
            << "Function: " << func << '\n';
  std::abort();
}

[[noreturn]] inline void error_message_impl(
  const std::string& message, const char* file, int line, const char* func) {
  std::cerr << "Error: " << message << '\n'
            << "File: " << file << '\n'
            << "Line: " << line << '\n'
            << "Function: " << func << '\n';
  std::abort();
}

}  // namespace detail

}  // namespace xcore

#ifdef DEBUG

// only enabled in DEBUG is defined
#define ASSERT(x)                                                             \
  do {                                                                        \
    if (!(x)) xcore::detail::assert_fail_impl(#x, __FILE__, __LINE__, __FUNCTION__); \
  } while (false)

#else

#define ASSERT(x) ((void)0)

#endif  // DEBUG


#define ENSURE(x)                                                             \
  do {                                                                        \
    if (!(x)) xcore::detail::assert_fail_impl(#x, __FILE__, __LINE__, __FUNCTION__); \
  } while (false)


#define ERROR_MESSAGE(stream_expr)                                             \
  do {                                                                         \
    std::ostringstream error_msg;                                              \
    error_msg << stream_expr;                                                  \
    xcore::detail::error_message_impl(error_msg.str(), __FILE__, __LINE__, __func__); \
  } while (0)

#endif  // XCORE_UTIL_ASSERT_H