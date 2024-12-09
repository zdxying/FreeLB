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

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "lbm/lattice_set.h"

namespace tempgenest {

struct codegenbase {

// often used chars
  static constexpr char _brace = '{';
  static constexpr char brace_ = '}';

  static constexpr char _bracket = '[';
  static constexpr char bracket_ = ']';

  static constexpr char _parenthesis = '(';
  static constexpr char parenthesis_ = ')';

  static constexpr char add = '+';
  static constexpr char sub = '-';
  static constexpr char mul = '*';
  static constexpr char div = '/';

  static constexpr char comma = ',';
  static constexpr char semicolon = ';';

  static constexpr char space = ' ';

  std::string _filename;

  std::vector<std::string> _args;
  std::vector<std::string> _argtypes;

  std::vector<std::string> _spectemplist;
  std::vector<unsigned int> _qlist;
};


struct rhoImplgen : public codegenbase {

  rhoImplgen(
  const std::string& filename,
  const std::vector<std::string>& spectemplist, 
  const std::vector<unsigned int>& qlist) {
    this->_filename = filename;
    this->_spectemplist = spectemplist;
    this->_qlist = qlist;
  }

  std::string _template = "template <typename T>";
  std::string _struct = "struct rhoImpl";
  std::string _function = "static void apply(const std::vector<T>& cell, T& rho_value)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename);
    for (unsigned int i = 0; i < this->_spectemplist.size(); ++i) {
      generate(file, this->_spectemplist[i], this->_qlist[i]);
    }
    file.close();
  }

  void generate(std::ofstream& file, std::string& spectemp, unsigned int q){
    file << _template << std::endl;
    file << _struct << "<" << spectemp << ">" << _brace << std::endl;
    file << _function << _brace << std::endl;

    // function body
    file << "rho_value = ";
    for (unsigned int i = 0; i < q; ++i) {
      file << "cell[" << i << "]";
      if (i < q - 1) file << add;
    }
    file << semicolon << std::endl;

    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


}  // namespace tempgenest

