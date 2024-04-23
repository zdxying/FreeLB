/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// ini_reader.h

#pragma once

#include <stdio.h>

#include <fstream>
#include <map>
#include <sstream>
#include <string>

// Read an INI file into easy-to-access name/value pairs.
// file example

/*
[int]            ; comment
intx = 25       ; comment
; comment
*/

/*
[section]
key=value
*/

// handle inline comments
// std::size_t end = line.find(';');
// if (end != std::string::npos) {
//   line = line.substr(0, end);
// }

class iniReader {
 private:
  std::map<std::string, std::map<std::string, std::string>> iniContent;
  // trim spaces from both ends of a string
  std::string trim(const std::string& str) {
    std::size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first) {
      return str;
    }
    std::size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
  }

 public:
  iniReader(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
      std::cerr << "Can't open file: " << filename << std::endl;
      exit(1);
    }
    std::string line, section;
    while (std::getline(file, line)) {
      // Ignore comments
      if (line[0] == ';') continue;

      if (line[0] == '[') {
        // read section
        section = line.substr(1, line.find(']') - 1);
      } else {
        // read key-value pair
        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '=')) {
          std::string value;
          if (std::getline(is_line, value)) {
            key = trim(key);
            value = trim(value);
            iniContent[section][key] = value;
          }
        }
      }
    }
  }

  // get value from section and key
  template <typename T>
  T getValue(const std::string& section, const std::string& name) {
    std::string value = iniContent[section][name];
    std::istringstream is_value(value);
    T val;
    is_value >> val;
    return val;
  }

  // get all values from a section
  template <typename T>
  void getVector(const std::string& section, std::vector<T>& vec) {
    for (auto& it : iniContent[section]) {
      std::istringstream is_value(it.second);
      T val;
      is_value >> val;
      vec.push_back(val);
    }
  }
};

// avoid unnecessary conversion
template <>
std::string iniReader::getValue<std::string>(const std::string& section,
                                             const std::string& name) {
  return iniContent[section][name];
}