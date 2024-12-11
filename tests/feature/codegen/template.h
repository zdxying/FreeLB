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

namespace tempgen {


struct D2Q5 : public Basic_Lattice_Set<2, 5> {
  static std::string name() { return "D2Q5"; }
  static std::string name(const std::string& s) { return "D2Q5" + s; }
  static constexpr Vector<int, 2> c[5] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}};
};;
struct D2Q9 : public Basic_Lattice_Set<2, 9> {
  static std::string name() { return "D2Q9"; }
  static std::string name(const std::string& s) { return "D2Q9" + s; }
  static constexpr Vector<int, 2> c[9] = {{0, 0}, {1, 0},   {-1, 0}, {0, 1}, {0, -1},
                                          {1, 1}, {-1, -1}, {1, -1}, {-1, 1}};
};
struct D3Q7 : public Basic_Lattice_Set<3, 7> {
  static std::string name() { return "D3Q7"; }
  static std::string name(const std::string& s) { return "D3Q7" + s; }
  static constexpr Vector<int, 3> c[7] = {{0, 0, 0},  {1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
                                          {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
};
struct D3Q15 : public Basic_Lattice_Set<3, 15> {
  static std::string name() { return "D3Q15"; }
  static std::string name(const std::string& s) { return "D3Q15" + s; }
  static constexpr Vector<int, 3> c[15] = {
    {0, 0, 0},   {1, 0, 0},  {-1, 0, 0},  {0, 1, 0},    {0, -1, 0},
    {0, 0, 1},   {0, 0, -1}, {1, 1, 1},   {-1, -1, -1}, {1, 1, -1},
    {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}, {-1, 1, 1},   {1, -1, -1}};
};
struct D3Q19 : public Basic_Lattice_Set<3, 19> {
  static std::string name() { return "D3Q19"; }
  static std::string name(const std::string& s) { return "D3Q19" + s; }
  static constexpr Vector<int, 3> c[19] = {
    {0, 0, 0},  {1, 0, 0},   {-1, 0, 0}, {0, 1, 0},   {0, -1, 0}, {0, 0, 1},   {0, 0, -1},
    {1, 1, 0},  {-1, -1, 0}, {1, 0, 1},  {-1, 0, -1}, {0, 1, 1},  {0, -1, -1}, {1, -1, 0},
    {-1, 1, 0}, {1, 0, -1},  {-1, 0, 1}, {0, 1, -1},  {0, -1, 1}};
};
struct D3Q27 : public Basic_Lattice_Set<3, 27> {
  static std::string name() { return "D3Q27"; }
  static std::string name(const std::string& s) { return "D3Q27" + s; }
  static constexpr Vector<int, 3> c[27] = {
    {0, 0, 0},  {1, 0, 0},    {-1, 0, 0}, {0, 1, 0},   {0, -1, 0}, {0, 0, 1},
    {0, 0, -1},  // 0-6
    {1, 1, 0},  {-1, -1, 0},  {1, 0, 1},  {-1, 0, -1}, {0, 1, 1},  {0, -1, -1},
    {1, -1, 0}, {-1, 1, 0},   {1, 0, -1}, {-1, 0, 1},  {0, 1, -1}, {0, -1, 1},  // 7-18
    {1, 1, 1},  {-1, -1, -1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1},
    {-1, 1, 1}, {1, -1, -1}}; 
};

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

  // cellname(LatSet::name("<T>"), "T", "TypePack")
  static std::string cellname(const std::string& LatSet, const std::string& T, const std::string& TypePack) {
    return "CELL<" + T + ", " + LatSet + ", " + TypePack + ">";
  }
};


// template <typename CELLTYPE, bool WriteToField>
// struct rhoImpl {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   using GenericRho = typename CELL::GenericRho;

//   __any__ static inline void apply(CELL& cell, T& rho_value) {
//     rho_value = T{};
//     for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
//     if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
//   }
// };
template<typename... LatSets>
struct rhoImplgen : public codegenbase {
  rhoImplgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack, bool WriteToField>";
  std::string _struct = "struct rhoImpl";

  std::string _function = "__any__ static inline void apply(CELLTYPE& cell, T& rho_value)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") + ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";" + 
                         "\nusing GenericRho = typename CELLTYPE::GenericRho;";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ", WriteToField>" << _brace << std::endl;
    file << _using << std::endl;
    file << std::endl;
    file << _function << _brace << std::endl;

    // function body
    file << "rho_value = ";
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      file << "cell[" << i << "]";
      if (i < LatSet::q - 1) file << add;
    }
    file << semicolon << std::endl;
    file << "if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;" << std::endl;

    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE, typename SOURCE, bool WriteToField>
// struct sourcerhoImpl {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   using GenericRho = typename CELL::GenericRho;

//   __any__ static inline void apply(CELL& cell, const T source, T& rho_value) {
//     rho_value = T{};
//     for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
//     // fOmega: avoid lattice artifact
//     rho_value += source * T{0.5};
//     if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
//   }
// };

template<typename... LatSets>
struct sourcerhoImplgen : public codegenbase {
  sourcerhoImplgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack, typename SOURCE, bool WriteToField>";
  std::string _struct = "struct sourcerhoImpl";

  std::string _function = "__any__ static inline void apply(CELLTYPE& cell, const T source, T& rho_value)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }
  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") + ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";" + 
                         "\nusing GenericRho = typename CELLTYPE::GenericRho;";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ", SOURCE, WriteToField>" << _brace << std::endl;
    file << _using << std::endl;
    file << std::endl;
    file << _function << _brace << std::endl;

    // function body
    file << "rho_value = ";
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      file << "cell[" << i << "]" << add;
    }
    file << "source * T{0.5}" << semicolon << std::endl;
    file << "if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;" << std::endl;

    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE, bool WriteToField>
// struct UImpl {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   using GenericRho = typename CELL::GenericRho;

//   __any__ static inline void apply(CELL& cell, Vector<T, LatSet::d>& u_value) {
//     u_value.clear();
//     T rho_value{};
//     for (unsigned int i = 0; i < LatSet::q; ++i) {
//       rho_value += cell[i];
//       u_value += latset::c<LatSet>(i) * cell[i];
//     }
//     u_value /= rho_value;
//     if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
//   }
// };
template<typename... LatSets>
struct UImplgen : public codegenbase {
  UImplgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack, bool WriteToField>";
  std::string _struct = "struct UImpl";

  std::string _function = "__any__ static inline void apply(CELLTYPE& cell, Vector<T, LatSet::d>& u_value)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") + ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";" + 
                         "\nusing GenericRho = typename CELLTYPE::GenericRho;";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ", WriteToField>" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _function << _brace << std::endl;

    // function body
    // rho
    file << "const T rho_value = ";
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      file << "cell[" << i << "]";
      if (i < LatSet::q - 1) file << add;
    }
    file << semicolon << std::endl;
    // u
    for (unsigned int dim = 0; dim < LatSet::d; ++dim){
      bool first = true;
      file << "u_value[" << dim <<"] = (";
      for (unsigned int i = 1; i < LatSet::q; ++i) {
        if (LatSet::c[i][dim] == 1) {
          if (first) {file << "cell[" << i << "]"; first = false;}
          else file << add << "cell[" << i << "]";
        } else if (LatSet::c[i][dim] == -1) {
          file << sub << "cell[" << i << "]"; first = false;
        }
      }
      file << ") / rho_value;" << std::endl;
    }
    
    file << "if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;" << std::endl;

    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE, typename ForceScheme, bool WriteToField, unsigned int dir>
// struct forceUImpl {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   using GenericRho = typename CELL::GenericRho;
//   static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;
//   __any__ static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value) {
//     u_value.clear();
//     T rho_value{};
//     for (unsigned int i = 0; i < LatSet::q; ++i) {
//       rho_value += cell[i];
//       u_value += latset::c<LatSet>(i) * cell[i];
//     }
//     u_value += f_alpha * T{0.5};
//     u_value /= rho_value;
//     if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
//   }
//   // for scalar force
//   __any__ static inline void apply(CELL& cell, const T f, Vector<T, LatSet::d>& u_value) {
//     u_value.clear();
//     T rho_value{};
//     for (unsigned int i = 0; i < LatSet::q; ++i) {
//       rho_value += cell[i];
//       u_value += latset::c<LatSet>(i) * cell[i];
//     }
//     u_value[scalardir] += f * T{0.5};
//     u_value /= rho_value;
//     if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
//   }
// };


template<typename... LatSets>
struct forceUImplgen : public codegenbase {
  forceUImplgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>";
  std::string _struct = "struct forceUImpl";

  std::string _staticconstexpr = "static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;";
  std::string _function0 = "__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value)";
  std::string _function1 = "__any__ static inline void apply(CELLTYPE& cell, const T f, Vector<T, LatSet::d>& u_value)";

  void generateAll(){
    std::ofstream file;
    // open and write from the end of the file
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") +  ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";" + 
                         "\nusing GenericRho = typename CELLTYPE::GenericRho;";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ", ForceScheme, WriteToField, dir>" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _staticconstexpr << "\n" << std::endl;

    // function 0
    file << _function0 << _brace << std::endl;
    // rho
    file << "const T rho_value = ";
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      file << "cell[" << i << "]";
      if (i < LatSet::q - 1) file << add;
    }
    file << semicolon << std::endl;
    // u
    for (unsigned int dim = 0; dim < LatSet::d; ++dim) {
      bool first = true;
      file << "u_value[" << dim <<"] = (";
      for (unsigned int i = 1; i < LatSet::q; ++i) {
        if (LatSet::c[i][dim] == 1) {
          if (first) {file << "cell[" << i << "]"; first = false;}
          else file << add << "cell[" << i << "]";
        } else if (LatSet::c[i][dim] == -1) {
          file << sub << "cell[" << i << "]"; first = false;
        }
      }
      file << " + f_alpha[" << dim << "]*T{0.5}) / rho_value;" << std::endl;
    }
    file << "if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;" << std::endl;
    // end of function0
    file << brace_ << std::endl;

    // function 1
    file << _function1 << _brace << std::endl;
    // rho
    file << "const T rho_value = ";
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      file << "cell[" << i << "]";
      if (i < LatSet::q - 1) file << add;
    }
    file << semicolon << std::endl;
    // u
    for (unsigned int dim = 0; dim < LatSet::d; ++dim) {
      bool first = true;
      file << "u_value[" << dim <<"] = ";
      for (unsigned int i = 1; i < LatSet::q; ++i) {
        if (LatSet::c[i][dim] == 1) {
          if (first) {file << "cell[" << i << "]"; first = false;}
          else file << add << "cell[" << i << "]";
        } else if (LatSet::c[i][dim] == -1) {
          file << sub << "cell[" << i << "]"; first = false;
        }
      }
      file << ";" << std::endl;
    }
    file << "u_value[scalardir] += f * T{0.5};" << std::endl;
    for (unsigned int dim = 0; dim < LatSet::d; ++dim) {
      file << "u_value[" << dim << "] /= rho_value;" << std::endl;
    }
    
    file << "if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;" << std::endl;

    // end of function1
    file << brace_ << std::endl;

    // end of struct
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE, bool WriteToField>
// struct rhoUImpl {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   using GenericRho = typename CELL::GenericRho;

//   __any__ static void apply(CELL& cell, T& rho_value, Vector<T, LatSet::d>& u_value) {
//     rho_value = T{};
//     u_value.clear();
//     for (unsigned int i = 0; i < LatSet::q; ++i) {
//       rho_value += cell[i];
//       u_value += latset::c<LatSet>(i) * cell[i];
//     }
//     u_value /= rho_value;
//     if constexpr (WriteToField) {
//       cell.template get<GenericRho>() = rho_value;
//       cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
//     }
//   }
// };
template<typename... LatSets>
struct rhoUImplgen : public codegenbase {
  rhoUImplgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack, bool WriteToField>";
  std::string _struct = "struct rhoUImpl";

  std::string _function = "__any__ static void apply(CELLTYPE& cell, T& rho_value, Vector<T, LatSet::d>& u_value)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") +  ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";" + 
                         "\nusing GenericRho = typename CELLTYPE::GenericRho;";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ", WriteToField>" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _function << _brace << std::endl;

    // function body
    // rho
    file << "rho_value = ";
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      file << "cell[" << i << "]";
      if (i < LatSet::q - 1) file << add;
    }
    file << semicolon << std::endl;
    // u
    for (unsigned int dim = 0; dim < LatSet::d; ++dim){
      bool first = true;
      file << "u_value[" << dim <<"] = (";
      for (unsigned int i = 1; i < LatSet::q; ++i) {
        if (LatSet::c[i][dim] == 1) {
          if (first) {file << "cell[" << i << "]"; first = false;}
          else file << add << "cell[" << i << "]";
        } else if (LatSet::c[i][dim] == -1) {
          file << sub << "cell[" << i << "]"; first = false;
        }
      }
      file << ") / rho_value;" << std::endl;
    }
    
    file << "if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}" << std::endl;

    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE, typename ForceScheme, bool WriteToField, unsigned int dir>
// struct forcerhoUImpl {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   using GenericRho = typename CELL::GenericRho;
//   static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;
//   __any__ static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value, Vector<T, LatSet::d>& u_value) {
//     rho_value = T{};
//     u_value.clear();
//     for (unsigned int i = 0; i < LatSet::q; ++i) {
//       rho_value += cell[i];
//       u_value += latset::c<LatSet>(i) * cell[i];
//     }
//     u_value += f_alpha * T{0.5};
//     u_value /= rho_value;
//     if constexpr (WriteToField) {
//       cell.template get<GenericRho>() = rho_value;
//       cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
//     }
//   }
//   // for scalar force
//   __any__ static inline void apply(CELL& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value) {
//     rho_value = T{};
//     u_value.clear();
//     for (unsigned int i = 0; i < LatSet::q; ++i) {
//       rho_value += cell[i];
//       u_value += latset::c<LatSet>(i) * cell[i];
//     }
//     u_value[scalardir] += f * T{0.5};
//     u_value /= rho_value;
//     if constexpr (WriteToField) {
//       cell.template get<GenericRho>() = rho_value;
//       cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
//     }
//   }
// };
template<typename... LatSets>
struct forcerhoUImplgen : public codegenbase {
  forcerhoUImplgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>";
  std::string _struct = "struct forcerhoUImpl";

  std::string _staticconstexpr = "static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;";
  std::string _function0 = "__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value, Vector<T, LatSet::d>& u_value)";
  std::string _function1 = "__any__ static inline void apply(CELLTYPE& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value)";

  void generateAll(){
    std::ofstream file;
    // open and write from the end of the file
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") +  ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";" + 
                         "\nusing GenericRho = typename CELLTYPE::GenericRho;";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ", ForceScheme, WriteToField, dir>" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _staticconstexpr << "\n" << std::endl;

    // function 0
    file << _function0 << _brace << std::endl;
    // rho
    file << "rho_value = ";
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      file << "cell[" << i << "]";
      if (i < LatSet::q - 1) file << add;
    }
    file << semicolon << std::endl;
    // u
    for (unsigned int dim = 0; dim < LatSet::d; ++dim) {
      bool first = true;
      file << "u_value[" << dim <<"] = (";
      for (unsigned int i = 1; i < LatSet::q; ++i) {
        if (LatSet::c[i][dim] == 1) {
          if (first) {file << "cell[" << i << "]"; first = false;}
          else file << add << "cell[" << i << "]";
        } else if (LatSet::c[i][dim] == -1) {
          file << sub << "cell[" << i << "]"; first = false;
        }
      }
      file << " + f_alpha[" << dim << "]*T{0.5}) / rho_value;" << std::endl;
    }
    file << "if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}" << std::endl;
    // end of function0
    file << brace_ << std::endl;

    // function 1
    file << _function1 << _brace << std::endl;
    // rho
    file << "rho_value = ";
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      file << "cell[" << i << "]";
      if (i < LatSet::q - 1) file << add;
    }
    file << semicolon << std::endl;
    // u
    for (unsigned int dim = 0; dim < LatSet::d; ++dim) {
      bool first = true;
      file << "u_value[" << dim <<"] = ";
      for (unsigned int i = 1; i < LatSet::q; ++i) {
        if (LatSet::c[i][dim] == 1) {
          if (first) {file << "cell[" << i << "]"; first = false;}
          else file << add << "cell[" << i << "]";
        } else if (LatSet::c[i][dim] == -1) {
          file << sub << "cell[" << i << "]"; first = false;
        }
      }
      file << ";" << std::endl;
    }
    file << "u_value[scalardir] += f * T{0.5};" << std::endl;
    for (unsigned int dim = 0; dim < LatSet::d; ++dim) {
      file << "u_value[" << dim << "] /= rho_value;" << std::endl;
    }
    
    file << "if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}" << std::endl;

    // end of function1
    file << brace_ << std::endl;

    // end of struct
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE>
// struct Pi_ab_neq {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   __any__ static inline void apply(CELL& cell, const T rho, const Vector<T, LatSet::d>& u, 
//   std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor) {
//     unsigned int i{};
//     for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
//       for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
//         T value{};
//         for (unsigned int k = 0; k < LatSet::q; ++k) {
//           value += latset::c<LatSet>(k)[alpha] * latset::c<LatSet>(k)[beta] * cell[k];
//         }
//         // remove the equilibrium part: PI_ab^eq
//         value -= rho * u[alpha] * u[beta];
//         if (alpha == beta) value -= rho * LatSet::cs2;
//         tensor[i] = value;
//         ++i;
//       }
//     }
//   }
// };
template<typename... LatSets>
struct Pi_ab_neqgen : public codegenbase {
  Pi_ab_neqgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack>";
  std::string _struct = "struct Pi_ab_neq";

  std::string _function = "__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") + ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ">" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _function << _brace << std::endl;

    // function body
    unsigned int i{};
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        file << "tensor[" << i << "] = ";
        bool first = true;
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          int c_ab = LatSet::c[k][alpha] * LatSet::c[k][beta];
          if (c_ab == 1) {
            if (first) {file << "cell[" << k << "]"; first = false;}
            else file << add << "cell[" << k << "]";
          } else if (c_ab == -1) {
            file << sub << "cell[" << k << "]"; first = false;
          }
        }
        // remove the equilibrium part: PI_ab^eq
        file << sub << " rho*(u[" << alpha << "]*u[" << beta << "]";
        // value -= rho * u[alpha] * u[beta]; if (alpha == beta) value -= rho * LatSet::cs2;
        // use ADD here to include LatSet::cs2 in ()
        if (alpha == beta) file << add << "LatSet::cs2";
        file << ");" << std::endl;
        ++i;
      }
    }
    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE>
// struct forcePi_ab_neq {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   __any__ static inline void apply(CELL& cell, const T rho, const Vector<T, LatSet::d>& u, const Vector<T, LatSet::d>& f_alpha, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor) {
//     unsigned int i{};
//     // remove force term in u
//     Vector<T, LatSet::d> unew = u - f_alpha * T{0.5}; //(T{0.5} / rho);
//     for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
//       for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
//         T value{};
//         const T force = T{0.5} * (f_alpha[alpha] * unew[beta] + f_alpha[beta] * unew[alpha]);  //* rho
//         for (unsigned int k = 0; k < LatSet::q; ++k) {
//           value += latset::c<LatSet>(k)[alpha] * latset::c<LatSet>(k)[beta] * cell[k];
//         }
//         value += force;
//         // remove the equilibrium part: PI_ab^eq
//         value -= rho * unew[alpha] * unew[beta];
//         if (alpha == beta) value -= rho * LatSet::cs2;
//         tensor[i] = value;
//         ++i;
//       }
//     }
//   }
// };
template<typename... LatSets>
struct forcePi_ab_neqgen : public codegenbase {
  forcePi_ab_neqgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack>";
  std::string _struct = "struct forcePi_ab_neq";
  std::string _function = "__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u_f, const Vector<T, LatSet::d>& f_alpha, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") + ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ">" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _function << _brace << std::endl;

    // function body
    // remove force term in u
    file << "const Vector<T, LatSet::d> u = u_f - f_alpha * T{0.5};" << std::endl;
    unsigned int i{};
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        file << "tensor[" << i << "] = ";
        bool first = true;
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          int c_ab = LatSet::c[k][alpha] * LatSet::c[k][beta];
          if (c_ab == 1) {
            if (first) {file << "cell[" << k << "]"; first = false;}
            else file << add << "cell[" << k << "]";
          } else if (c_ab == -1) {
            file << sub << "cell[" << k << "]"; first = false;
          }
        }
        // add force term
        if (!first) file << add << space; 
        file << "T{0.5}*(f_alpha[" << alpha << "]*u[" << beta << "]+f_alpha[" << beta << "]*u[" << alpha << "]) ";
        // remove the equilibrium part: PI_ab^eq
        file << sub << " rho*(u[" << alpha << "]*u[" << beta << "]";
        // value -= rho * unew[alpha] * unew[beta]; if (alpha == beta) value -= rho * LatSet::cs2;
        // use ADD here to include LatSet::cs2 in ()
        if (alpha == beta) file << add << "LatSet::cs2";
        file << ");" << std::endl;
        ++i;
      }
    }
    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE>
// struct stress {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   __any__ static inline void apply(CELL& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& stress_tensor) {
//     unsigned int i{};
//     const T coeff = T{0.5} * cell.getOmega() - T{1};
//     for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
//       for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
//         T value{};
//         for (unsigned int k = 0; k < LatSet::q; ++k) {
//           value += latset::c<LatSet>(k)[alpha] * latset::c<LatSet>(k)[beta] * cell[k];
//         }
//         // remove the equilibrium part: PI_ab^eq
//         value -= rho * u[alpha] * u[beta];
//         if (alpha == beta) value -= rho * LatSet::cs2;
//         // multiply by the coefficient
//         value *= coeff;
//         stress_tensor[i] = value;
//         ++i;
//       }
//     }
//   }
// };
template<typename... LatSets>
struct stressgen : public codegenbase {
  stressgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack>";
  std::string _struct = "struct stress";

  std::string _function = "__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") + ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ">" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _function << _brace << std::endl;

    // function body
    // coeffient
    file << "const T coeff = T{0.5} * cell.getOmega() - T{1};" << std::endl;
    unsigned int i{};
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        file << "tensor[" << i << "] = (";
        bool first = true;
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          int c_ab = LatSet::c[k][alpha] * LatSet::c[k][beta];
          if (c_ab == 1) {
            if (first) {file << "cell[" << k << "]"; first = false;}
            else file << add << "cell[" << k << "]";
          } else if (c_ab == -1) {
            file << sub << "cell[" << k << "]"; first = false;
          }
        }
        // remove the equilibrium part: PI_ab^eq
        file << sub << "rho*(u[" << alpha << "]*u[" << beta << "]";
        // value -= rho * u[alpha] * u[beta]; if (alpha == beta) value -= rho * LatSet::cs2;
        // use ADD here to include LatSet::cs2 in ()
        if (alpha == beta) file << add << "LatSet::cs2";
        file << ")) * coeff;" << std::endl;
        ++i;
      }
    }
    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE>
// struct strainRate {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   __any__ static inline void apply(CELL& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& strain_rate_tensor) {
//     unsigned int i{};
//     const T coeff = T{-1.5} * cell.getOmega() / cell.template get<typename CELL::GenericRho>();
//     for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
//       for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
//         T value{};
//         for (unsigned int k = 0; k < LatSet::q; ++k) {
//           value += latset::c<LatSet>(k)[alpha] * latset::c<LatSet>(k)[beta] * cell[k];
//         }
//         // remove the equilibrium part: PI_ab^eq
//         value -= rho * u[alpha] * u[beta];
//         if (alpha == beta) value -= rho * LatSet::cs2;
//         // multiply by the coefficient
//         value *= coeff;
//         strain_rate_tensor[i] = value;
//         ++i;
//       }
//     }
//   }
// };
template<typename... LatSets>
struct strainRategen : public codegenbase {
  strainRategen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack>";
  std::string _struct = "struct strainRate";

  std::string _function = "__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") + ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ">" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _function << _brace << std::endl;

    // function body
    // coeffient
    file << "const T coeff = T{-1.5} * cell.getOmega() / cell.template get<typename CELLTYPE::GenericRho>();" << std::endl;
    unsigned int i{};
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        file << "tensor[" << i << "] = (";
        bool first = true;
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          int c_ab = LatSet::c[k][alpha] * LatSet::c[k][beta];
          if (c_ab == 1) {
            if (first) {file << "cell[" << k << "]"; first = false;}
            else file << add << "cell[" << k << "]";
          } else if (c_ab == -1) {
            file << sub << "cell[" << k << "]"; first = false;
          }
        }
        // remove the equilibrium part: PI_ab^eq
        file << sub << "rho*(u[" << alpha << "]*u[" << beta << "]";
        // value -= rho * u[alpha] * u[beta]; if (alpha == beta) value -= rho * LatSet::cs2;
        // use ADD here to include LatSet::cs2 in ()
        if (alpha == beta) file << add << "LatSet::cs2";
        file << ")) * coeff;" << std::endl;
        ++i;
      }
    }
    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELLTYPE, bool WriteToField>
// struct shearRateMagImpl {
//   using CELL = CELLTYPE;
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   __any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& strain_rate_tensor) {
//     T value{};
//     unsigned int i{};
//     for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
//       for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
//         T sq = strain_rate_tensor[i] * strain_rate_tensor[i];
//         if (alpha != beta) sq *= T{2};
//         value += sq;
//         ++i;
//       }
//     }
//     return std::sqrt(T{2} * value);
//   }
// };
template<typename... LatSets>
struct shearRateMagImplgen : public codegenbase {
  shearRateMagImplgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack>";
  std::string _struct = "struct shearRateMagImpl";

  std::string _function = "__any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& srtensor)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") + ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ">" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _function << _brace << std::endl;

    // function body
    file << "T value = ";
    unsigned int i{};
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        file << "srtensor[" << i << "]*srtensor[" << i << "]";
        if (alpha != beta) file << mul << "2";
        // remove the equilibrium part: PI_ab^eq
        file << space;
        if (i < LatSet::d*(LatSet::d+1)/2 - 1) file << add << space;
        ++i;
      }
    }
    file << semicolon << std::endl;
    file << "return std::sqrt(2 * value);" << std::endl;
    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};


// template <typename CELL>
// struct SecondOrderImpl {
//   using T = typename CELL::FloatType;
//   using LatSet = typename CELL::LatticeSet;
//   using CELLTYPE = CELL;
//   using GenericRho = typename CELL::GenericRho;
//   __any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u) {
//     const T u2 = u.getnorm2();
//     for (unsigned int k = 0; k < LatSet::q; ++k) {
//       const T uc = u * latset::c<LatSet>(k);
//       feq[k] = latset::w<LatSet>(k) * rho *
//       (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 - LatSet::InvCs2 * u2 * T{0.5});
//     }
//   }
// };
template<typename... LatSets>
struct SecondOrderImplgen : public codegenbase {
  SecondOrderImplgen(const std::string& filename) { _filename = filename; }

  std::string _template = "template <typename T, typename TypePack>";
  std::string _struct = "struct SecondOrderImpl";
  std::string _function = "__any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u)";

  void generateAll(){
    std::ofstream file;
    file.open(_filename, std::ios::out | std::ios::app);
    file << "\n//------------------------------------\n" << std::endl;
    (generate<LatSets>(file), ...);
    file << "\n//------------------------------------\n" << std::endl;
    file.close();
    std::cout << "Generated " << _filename << std::endl;
  }

  template<typename LatSet>
  void generate(std::ofstream& file){
    std::string _using = "using CELLTYPE = " + cellname(LatSet::name("<T>"), "T", "TypePack") + ";" +
                         "\nusing LatSet = " + LatSet::name("<T>") + ";";

    file << _template << std::endl;
    file << _struct << "<"<< cellname(LatSet::name("<T>"), "T", "TypePack") << ">" << _brace << std::endl;
    file << _using << "\n"<< std::endl;
    file << _function << _brace << std::endl;

    // function body

    file << "constexpr T InvCs2x = T{0.5} * LatSet::InvCs2;" << std::endl;
    // var0 = T{1} - LatSet::InvCs2 * T{0.5} * u2 
    file << "const T var0 = T{1} - (u[0]*u[0] + u[1]*u[1]";
    if (LatSet::d == 3) file << " + u[2]*u[2]";
    file << ") * InvCs2x;" << std::endl;

    // rho * w[k]
    file << "const T rhowk0 = latset::w<LatSet>(0) * rho;" << std::endl;
    std::vector<unsigned int> zerosvec(LatSet::q, 0);
    unsigned int zerosbuffer = 3;
    for (unsigned int i = 1; i < LatSet::q; ++i) {
      // find number of zeros in LatSet::c
      unsigned int zeros{};
      for (unsigned int dim = 0; dim < LatSet::d; ++dim) {
        if (LatSet::c[i][dim] == 0) zeros++;
      }
      zerosvec[i] = zerosvec[i-1];
      if (zeros != zerosbuffer) {
        file << "const T rhowk" << i << " = latset::w<LatSet>(" << i << ") * rho;" << std::endl;
        zerosbuffer = zeros;
        zerosvec[i] = i;
      }
    }

    // uc terms
    for (unsigned int i = 1; i < LatSet::q; i += 2) {
      // LatSet::InvCs2 * uc
      file << "const T InvCs2uck" << i << " = LatSet::InvCs2 * (";
      bool first = true;
      for (unsigned int dim = 0; dim < LatSet::d; ++dim) {
        if (LatSet::c[i][dim] == 1) {
          if (first) {file << "u[" << dim << "]"; first = false;}
          else file << add << "u[" << dim << "]";
        } else if (LatSet::c[i][dim] == -1) {
          file << sub << "u[" << dim << "]"; first = false;
        }
      }
      file << ");" << std::endl;
      // uc * uc * T{0.5} * LatSet::InvCs4
      file << "const T var0InvCs4uc2k" << i << " = var0 + InvCs2uck" << i << "*InvCs2uck" << i << "*T{0.5};" << std::endl;
    }

    file << "feq[0] = rhowk0 * var0;" << std::endl;
    for (unsigned int i = 1; i < LatSet::q; i += 2) {
      file << "feq[" << i   << "] = rhowk" << zerosvec[i] << " * (var0InvCs4uc2k" << i <<" + InvCs2uck" << i << ");" << std::endl;
      file << "feq[" << i+1 << "] = rhowk" << zerosvec[i] << " * (var0InvCs4uc2k" << i <<" - InvCs2uck" << i << ");" << std::endl;
    }

    // end of function and struct
    file << brace_ << std::endl;
    file << "};" <<"\n"<< std::endl;
  }
};
}  // namespace tempgen

