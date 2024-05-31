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

// vtkWriter.h

#pragma once

#include <any>

#include "data_struct/Vector.h"
#include "data_struct/lattice.h"
#include "io/basic_writer.h"

// vtk format:
// dataType is one of the types: bit, unsigned_char, char, unsigned_short,
// short, unsigned_int, int, unsigned_long, long, float, or double
class AbstractFieldWriter {
 public:
  virtual void write(std::ofstream &f) = 0;
};

namespace vtkWriter {

template <typename ArrayType>
class FlagWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const ArrayType &Field;
  std::size_t Size;

 public:
  FlagWriter(std::string name, const ArrayType &f)
      : varname(name), Field(f), Size(Field.size()) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    ss << "SCALARS " << varname << " int" << std::endl;
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Size; ++i) {
      ss << static_cast<int>(Field[i]) << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};
template <typename FIELDTYPE>
class FieldFlagWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const FIELDTYPE *Field;
  const int Size;

 public:
  FieldFlagWriter(std::string name, const FIELDTYPE *f, int size)
      : varname(name), Field(f), Size(size) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    ss << "SCALARS " << varname << " int" << std::endl;
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Size; ++i) {
      ss << static_cast<int>(Field[i]) << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};

template <template <typename> class ArrayType, typename T>
class ScalarWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const ArrayType<T> &Field;
  std::size_t Size;

 public:
  ScalarWriter(std::string name, const ArrayType<T> &f)
      : varname(name), Field(f), Size(Field.size()) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<T, double>::value) {
      ss << "SCALARS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<T, float>::value) {
      ss << "SCALARS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<T, int>::value) {
      ss << "SCALARS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: ScalarWriter: unsupported type" << std::endl;
      exit(1);
    }
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Size; ++i) {
      ss << Field[i] << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};
template <typename FIELDTYPE>
class FieldScalarWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const FIELDTYPE *Field;
  const int Size;

 public:
  FieldScalarWriter(std::string name, const FIELDTYPE *f, int size)
      : varname(name), Field(f), Size(size) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<FIELDTYPE, double>::value) {
      ss << "SCALARS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, float>::value) {
      ss << "SCALARS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, int>::value) {
      ss << "SCALARS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: FieldScalarWriter: unsupported type" << std::endl;
      exit(1);
    }
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Size; ++i) {
      ss << Field[i] << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};

template <template <typename> class ArrayType, typename T>
class PhysScalarWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const ArrayType<T> &Field;
  std::size_t Size;
  const AbstractConverter<T> &Conv;

 public:
  PhysScalarWriter(std::string name, const ArrayType<T> &f,
                   const AbstractConverter<T> &conv)
      : varname(name), Field(f), Size(f.size()), Conv(conv) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<T, double>::value) {
      ss << "SCALARS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<T, float>::value) {
      ss << "SCALARS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<T, int>::value) {
      ss << "SCALARS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: PhysScalarWriter: unsupported type" << std::endl;
      exit(1);
    }
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Size; ++i) {
      ss << Conv.getPhysRho(Field[i]) << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};
template <typename FIELDTYPE>
class PhysFieldScalarWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const FIELDTYPE *Field;
  const int Size;
  const AbstractConverter<FIELDTYPE> &Conv;

 public:
  PhysFieldScalarWriter(std::string name, const FIELDTYPE *f, int size,
                        const AbstractConverter<FIELDTYPE> &conv)
      : varname(name), Field(f), Size(size), Conv(conv) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<FIELDTYPE, double>::value) {
      ss << "SCALARS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, float>::value) {
      ss << "SCALARS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, int>::value) {
      ss << "SCALARS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: PhysFieldScalarWriter: unsupported type" << std::endl;
      exit(1);
    }
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Size; ++i) {
      ss << Conv.getPhysRho(Field[i]) << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};

template <typename FIELDTYPE, unsigned int D>
class FieldVectorWriter_SOA : public AbstractFieldWriter {
 private:
  std::string varname;
  std::array<FIELDTYPE *, D> &Field;
  const int Size;

 public:
  template <typename... Args>
  FieldVectorWriter_SOA(std::string name, int size, Args... args)
      : varname(name), Size(size), Field{args...} {}
  FieldVectorWriter_SOA(std::string name, int size, std::array<FIELDTYPE *, D> &f)
      : varname(name), Size(size), Field(f) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<FIELDTYPE, double>::value) {
      ss << "VECTORS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, float>::value) {
      ss << "VECTORS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, int>::value) {
      ss << "VECTORS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: FieldVectorWriter_SOA: unsupported type" << std::endl;
      exit(1);
    }
    for (int i = 0; i < Size; ++i) {
      if constexpr (D == 2) {
        ss << Field[0][i] << " " << Field[1][i] << " " << 0 << " ";
      } else if constexpr (D == 3) {
        ss << Field[0][i] << " " << Field[1][i] << " " << Field[2][i] << " ";
      }
    }
    ss << std::endl;
    f << ss.str();
  }
};

template <typename FIELDTYPE, unsigned int D>
class FieldVectorWriter_AOS : public AbstractFieldWriter {
 private:
  std::string varname;
  const Vector<FIELDTYPE, D> *Field;
  const int Size;

 public:
  template <typename... Args>
  FieldVectorWriter_AOS(std::string name, const Vector<FIELDTYPE, D> *f, int size)
      : varname(name), Field(f), Size(size) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<FIELDTYPE, double>::value) {
      ss << "VECTORS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, float>::value) {
      ss << "VECTORS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, int>::value) {
      ss << "VECTORS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: FieldVectorWriter_AOS: unsupported type" << std::endl;
      exit(1);
    }
    for (int i = 0; i < Size; ++i) {
      if constexpr (D == 2) {
        ss << Field[i][0] << " " << Field[i][1] << " " << 0 << " ";
      } else if constexpr (D == 3) {
        ss << Field[i][0] << " " << Field[i][1] << " " << Field[i][2] << " ";
      }
    }
    ss << std::endl;
    f << ss.str();
  }
};

template <typename T, unsigned int D>
class PhysVelocityWriter_AOS : public AbstractFieldWriter {
 private:
  std::string varname;
  const GenericArray<Vector<T, D>> &Field;
  std::size_t Size;
  const AbstractConverter<T> &Conv;

 public:
  template <typename... Args>
  PhysVelocityWriter_AOS(std::string name, const GenericArray<Vector<T, D>> &f,
                         const AbstractConverter<T> &conv)
      : varname(name), Field(f), Size(Field.size()), Conv(conv) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<T, double>::value) {
      ss << "VECTORS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<T, float>::value) {
      ss << "VECTORS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<T, int>::value) {
      ss << "VECTORS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: FieldVectorWriter_AOS: unsupported type" << std::endl;
      exit(1);
    }
    for (int i = 0; i < Size; ++i) {
      if constexpr (D == 2) {
        ss << Conv.getPhysU(Field[i][0]) << " " << Conv.getPhysU(Field[i][1]) << " " << 0
           << " ";
      } else if constexpr (D == 3) {
        ss << Conv.getPhysU(Field[i][0]) << " " << Conv.getPhysU(Field[i][1]) << " "
           << Conv.getPhysU(Field[i][2]) << " ";
      }
    }
    ss << std::endl;
    f << ss.str();
  }
};
template <typename FIELDTYPE, unsigned int D>
class PhysVelocityFieldWriter_AOS : public AbstractFieldWriter {
 private:
  std::string varname;
  const Vector<FIELDTYPE, D> *Field;
  const int Size;
  const AbstractConverter<FIELDTYPE> &Conv;

 public:
  template <typename... Args>
  PhysVelocityFieldWriter_AOS(std::string name, const Vector<FIELDTYPE, D> *f, int size,
                              const AbstractConverter<FIELDTYPE> &conv)
      : varname(name), Field(f), Size(size), Conv(conv) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<FIELDTYPE, double>::value) {
      ss << "VECTORS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, float>::value) {
      ss << "VECTORS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, int>::value) {
      ss << "VECTORS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: FieldVectorWriter_AOS: unsupported type" << std::endl;
      exit(1);
    }
    for (int i = 0; i < Size; ++i) {
      if constexpr (D == 2) {
        ss << Conv.getPhysU(Field[i][0]) << " " << Conv.getPhysU(Field[i][1]) << " " << 0
           << " ";
      } else if constexpr (D == 3) {
        ss << Conv.getPhysU(Field[i][0]) << " " << Conv.getPhysU(Field[i][1]) << " "
           << Conv.getPhysU(Field[i][2]) << " ";
      }
    }
    ss << std::endl;
    f << ss.str();
  }
};

// ---------------------------------
// for unstructured grid
template <typename T>
class UnStruFlagWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const std::vector<T> &Field;

 public:
  UnStruFlagWriter(std::string name, std::vector<T> &f) : varname(name), Field(f) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    ss << "SCALARS " << varname << " int" << std::endl;
    ss << "LOOKUP_TABLE default" << std::endl;
    for (auto field : Field) {
      ss << static_cast<int>(field) << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};

template <typename T>
class UnStruScalarWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const std::vector<T> &Field;

 public:
  UnStruScalarWriter(std::string name, std::vector<T> &f) : varname(name), Field(f) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<T, double>::value) {
      ss << "SCALARS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<T, float>::value) {
      ss << "SCALARS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<T, int>::value) {
      ss << "SCALARS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: ScalarWriter: unsupported type" << std::endl;
      exit(1);
    }
    ss << "LOOKUP_TABLE default" << std::endl;
    for (auto field : Field) {
      ss << field << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};


}  // namespace vtkWriter

template <typename FIELDTYPE>
class FieldFlagWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const std::vector<FIELDTYPE> &field;
  const std::vector<int> &index;

 public:
  FieldFlagWriter(std::string name, std::vector<FIELDTYPE> &f, std::vector<int> &idx)
      : varname(name), field(f), index(idx) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    ss << "SCALARS " << varname << " int" << std::endl;
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < index.size(); ++i) {
      if (index[i] != -1)
        ss << static_cast<int>(field[index[i]]) << " ";
      else
        ss << -2 << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};

template <typename FIELDTYPE>
class FieldScalarWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const std::vector<FIELDTYPE> &field;
  const std::vector<int> &index;

 public:
  FieldScalarWriter(std::string name, std::vector<FIELDTYPE> &f, std::vector<int> &idx)
      : varname(name), field(f), index(idx) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<FIELDTYPE, double>::value) {
      ss << "SCALARS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, float>::value) {
      ss << "SCALARS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, int>::value) {
      ss << "SCALARS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: FieldScalarWriter: unsupported type" << std::endl;
      exit(1);
    }
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < index.size(); ++i) {
      if (index[i] != -1)
        ss << field[index[i]] << " ";
      else
        ss << 0 << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
};

template <typename FIELDTYPE, unsigned int D>
class FieldVectorWriter : public AbstractFieldWriter {
 private:
  std::string varname;
  const std::vector<Vector<FIELDTYPE, D>> &field;
  const std::vector<int> &index;

 public:
  FieldVectorWriter(std::string name, std::vector<Vector<FIELDTYPE, D>> &f,
                    std::vector<int> &idx)
      : varname(name), field(f), index(idx) {}
  void write(std::ofstream &f) override {
    std::stringstream ss;
    if constexpr (std::is_same<FIELDTYPE, double>::value) {
      ss << "VECTORS " << varname << " double" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, float>::value) {
      ss << "VECTORS " << varname << " float" << std::endl;
    } else if constexpr (std::is_same<FIELDTYPE, int>::value) {
      ss << "VECTORS " << varname << " int" << std::endl;
    } else {
      std::cout << "ERROR: FieldVectorWriter: unsupported type" << std::endl;
      exit(1);
    }
    for (int i = 0; i < index.size(); ++i) {
      if (index[i] != -1) {
        if constexpr (D == 2) {
          ss << field[index[i]][0] << " " << field[index[i]][1] << " " << 0 << " ";
        } else {
          ss << field[index[i]][0] << " " << field[index[i]][1] << " "
             << field[index[i]][2] << " ";
        }
      } else {
        ss << 0 << " " << 0 << " " << 0 << " ";
      }
    }
    ss << std::endl;
    f << ss.str();
  }
};

template <typename T, unsigned int D>
class vtkStruPointsWriter {
 private:
  int _Nx;
  int _Ny;
  int _Nz;
  T VoxelSize;
  Vector<T, D> _Min;
  std::string _dirname = "./vtkoutput";
  std::string _filename;
  // field data to write
  std::vector<AbstractFieldWriter *> _FieldWriters;

 public:
  vtkStruPointsWriter(std::string filename, Geometry<T, D> &geo)
      : _filename(filename), VoxelSize(geo.getVoxelSize()), _Min(geo.getMin()),
        _Nx(geo.getNx()), _Ny(geo.getNy()), _Nz(geo.getNz()) {}
  vtkStruPointsWriter(std::string filename, T voxsize, Vector<T, D> Min, int Nx, int Ny,
                      int Nz = 1)
      : _filename(filename), VoxelSize(voxsize), _Min(Min), _Nx(Nx), _Ny(Ny), _Nz(Nz) {}

  void addtoWriteList(AbstractFieldWriter *writer) { _FieldWriters.push_back(writer); }
  template <typename... Args>
  void addtoWriteList(AbstractFieldWriter *writer, Args... args) {
    _FieldWriters.push_back(writer);
    addtoWriteList(args...);
  }
  void setDirname(std::string dirname) { _dirname = dirname; }

  void Write() {
    DirCreator::Create_Dir(_dirname);
    std::string fullName = "./vtkoutput/" + _filename + ".vtk";
    std::ofstream f(fullName.c_str());
    writeHeader(f);
    // write
    for (auto writer : _FieldWriters) {
      writer->write(f);
    }
    f.close();
  }
  void Write(int step) {
    DirCreator::Create_Dir(_dirname);
    std::string fullName = "./vtkoutput/" + _filename + std::to_string(step) + ".vtk";
    std::ofstream f(fullName.c_str());
    writeHeader(f);
    // write
    for (auto writer : _FieldWriters) {
      writer->write(f);
    }
    f.close();
  }

  void writeHeader(std::ofstream &f) {
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "Voxels" << std::endl;
    f << "ASCII" << std::endl;

    f << "DATASET STRUCTURED_POINTS" << std::endl;
    f << "DIMENSIONS " << _Nx << " " << _Ny << " " << _Nz << std::endl;
    if constexpr (D == 2) {
      f << "ORIGIN " << _Min[0] + VoxelSize * T(0.5) << " "
        << _Min[1] + VoxelSize * T(0.5) << " 0" << std::endl;
      f << "SPACING " << VoxelSize << " " << VoxelSize << " 1" << std::endl;
    } else {
      f << "ORIGIN " << _Min[0] + VoxelSize * T(0.5) << " "
        << _Min[1] + VoxelSize * T(0.5) << " " << _Min[2] + VoxelSize * T(0.5)
        << std::endl;
      f << "SPACING " << VoxelSize << " " << VoxelSize << " " << VoxelSize << std::endl;
    }
    f << "POINT_DATA " << _Nx * _Ny * _Nz << std::endl;
  }
};

template <typename T, unsigned int D>
class vtkUnStruGridWriter {
 private:
  std::string _dirname = "./vtkoutput";
  std::string _filename;
  // field data to write
  std::vector<AbstractFieldWriter *> _FieldWriters;

 public:
  vtkUnStruGridWriter(std::string filename) : _filename(filename) {}

  void addtoWriteList(AbstractFieldWriter *writer) { _FieldWriters.push_back(writer); }
  template <typename... Args>
  void addtoWriteList(AbstractFieldWriter *writer, Args... args) {
    _FieldWriters.push_back(writer);
    addtoWriteList(args...);
  }
  void setDirname(std::string dirname) { _dirname = dirname; }

  void Write(std::vector<BasicTree<T, D> *> &leafs) {
    DirCreator::Create_Dir(_dirname);
    std::string fullName = "./vtkoutput/" + _filename + ".vtk";
    std::ofstream f(fullName.c_str());
    writeHeader(f);
    writeGeometry(f, leafs);
    // write
    for (auto writer : _FieldWriters) {
      writer->write(f);
    }
    f.close();
  }
  void Write(std::vector<BasicTree<T, D> *> &leafs, int step) {
    DirCreator::Create_Dir(_dirname);
    std::string fullName = "./vtkoutput/" + _filename + std::to_string(step) + ".vtk";
    std::ofstream f(fullName.c_str());
    writeHeader(f);
    writeGeometry(f, leafs);
    // write
    for (auto writer : _FieldWriters) {
      writer->write(f);
    }
    f.close();
  }

  void writeHeader(std::ofstream &f) {
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "tree" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET UNSTRUCTURED_GRID" << std::endl;
  }

  void writeGeometry(std::ofstream &f, std::vector<BasicTree<T, D> *> &leafs) {
    std::stringstream points;
    std::stringstream cells;
    std::stringstream cell_types;
    std::stringstream celldataheader;

    if constexpr (D == 2) {
      if constexpr (std::is_same<T, double>::value) {
        points << "POINTS " << leafs.size() * 4 << " double" << std::endl;
      } else if constexpr (std::is_same<T, float>::value) {
        points << "POINTS " << leafs.size() * 4 << " float" << std::endl;
      } else if constexpr (std::is_same<T, int>::value) {
        points << "POINTS " << leafs.size() * 4 << " int" << std::endl;
      } else {
        std::cout << "ERROR: vtkUnStruGridWriter: unsupported type" << std::endl;
        exit(1);
      }
    } else if constexpr (D == 3) {
      if constexpr (std::is_same<T, double>::value) {
        points << "POINTS " << leafs.size() * 8 << " double" << std::endl;
      } else if constexpr (std::is_same<T, float>::value) {
        points << "POINTS " << leafs.size() * 8 << " float" << std::endl;
      } else if constexpr (std::is_same<T, int>::value) {
        points << "POINTS " << leafs.size() * 8 << " int" << std::endl;
      } else {
        std::cout << "ERROR: vtkUnStruGridWriter: unsupported type" << std::endl;
        exit(1);
      }
    }
    if constexpr (D == 2) {
      cells << "CELLS " << leafs.size() << " " << leafs.size() * 5 << std::endl;
    } else if constexpr (D == 3) {
      cells << "CELLS " << leafs.size() << " " << leafs.size() * 9 << std::endl;
    }
    cell_types << "CELL_TYPES " << leafs.size() << std::endl;
    T rad;
    int i = 0;
    if constexpr (D == 2) {
      Vector<T, 2> center;
      Vector<T, 2> pt;
      for (auto leaf : leafs) {
        center = leaf->getCenter();
        rad = leaf->getRadius();

        pt = center + Vector<T, 2>(-rad, -rad);
        points << pt[0] << " " << pt[1] << " 0 ";
        pt = center + Vector<T, 2>(rad, -rad);
        points << pt[0] << " " << pt[1] << " 0 ";
        pt = center + Vector<T, 2>(-rad, rad);
        points << pt[0] << " " << pt[1] << " 0 ";
        pt = center + Vector<T, 2>(rad, rad);
        points << pt[0] << " " << pt[1] << " 0 \n";

        cells << "4 " << i << " " << i + 1 << " " << i + 2 << " " << i + 3 << " \n";
        i += 4;
        cell_types << 8 << "\n";
      }
    } else if constexpr (D == 3) {
      Vector<T, 3> center;
      Vector<T, 3> pt;
      for (auto leaf : leafs) {
        center = leaf->getCenter();
        rad = leaf->getRadius();

        pt = center + Vector<T, 3>(-rad, -rad, -rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(rad, -rad, -rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(rad, rad, -rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(-rad, rad, -rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(-rad, -rad, rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(rad, -rad, rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(rad, rad, rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(-rad, rad, rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " \n";

        cells << "8 " << i << " " << i + 1 << " " << i + 2 << " " << i + 3 << " " << i + 4
              << " " << i + 5 << " " << i + 6 << " " << i + 7 << " \n";
        i += 8;
        cell_types << 11 << "\n";
      }
    }
    celldataheader << "CELL_DATA " << leafs.size() << std::endl;
    f << points.str() << cells.str() << cell_types.str() << celldataheader.str();
  }
};

template <typename T, unsigned int D>
class vtkPolyWriter {
 private:
  std::string _dirname = "./vtkoutput";
  std::string _filename;
  // field data to write
  std::vector<AbstractFieldWriter *> _FieldWriters;

 public:
  vtkPolyWriter(std::string filename) : _filename(filename) {}

  void addtoWriteList(AbstractFieldWriter *writer) { _FieldWriters.push_back(writer); }
  template <typename... Args>
  void addtoWriteList(AbstractFieldWriter *writer, Args... args) {
    _FieldWriters.push_back(writer);
    addtoWriteList(args...);
  }
  void setDirname(std::string dirname) { _dirname = dirname; }

  void Write(std::vector<BasicTree<T, D> *> &leafs) {
    DirCreator::Create_Dir(_dirname);
    std::string fullName = "./vtkoutput/" + _filename + ".vtk";
    std::ofstream f(fullName.c_str());
    writeHeader(f);
    writeGeometry(f, leafs);
    // write
    for (auto writer : _FieldWriters) {
      writer->write(f);
    }
    f.close();
  }
  void Write(std::vector<BasicTree<T, D> *> &leafs, int step) {
    DirCreator::Create_Dir(_dirname);
    std::string fullName = "./vtkoutput/" + _filename + std::to_string(step) + ".vtk";
    std::ofstream f(fullName.c_str());
    writeHeader(f);
    writeGeometry(f, leafs);
    // write
    for (auto writer : _FieldWriters) {
      writer->write(f);
    }
    f.close();
  }

  void writeHeader(std::ofstream &f) {
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "polydata" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET POLYDATA" << std::endl;
  }

  void writeGeometry(std::ofstream &f, std::vector<BasicTree<T, D> *> &leafs) {
    std::stringstream points;
    std::stringstream cells;
    std::stringstream celldataheader;

    if constexpr (D == 2) {
      if constexpr (std::is_same<T, double>::value) {
        points << "POINTS " << leafs.size() * 4 << " double" << std::endl;
      } else if constexpr (std::is_same<T, float>::value) {
        points << "POINTS " << leafs.size() * 4 << " float" << std::endl;
      } else if constexpr (std::is_same<T, int>::value) {
        points << "POINTS " << leafs.size() * 4 << " int" << std::endl;
      } else {
        std::cout << "ERROR: vtkPolyWriter: unsupported type" << std::endl;
        exit(1);
      }
    } else if constexpr (D == 3) {
      if constexpr (std::is_same<T, double>::value) {
        points << "POINTS " << leafs.size() * 8 << " double" << std::endl;
      } else if constexpr (std::is_same<T, float>::value) {
        points << "POINTS " << leafs.size() * 8 << " float" << std::endl;
      } else if constexpr (std::is_same<T, int>::value) {
        points << "POINTS " << leafs.size() * 8 << " int" << std::endl;
      } else {
        std::cout << "ERROR: vtkPolyWriter: unsupported type" << std::endl;
        exit(1);
      }
    }
    if constexpr (D == 2) {
      cells << "POLYGONS " << leafs.size() << " " << leafs.size() * 5 << std::endl;
    } else if constexpr (D == 3) {
      cells << "POLYGONS " << leafs.size() << " " << leafs.size() * 9 << std::endl;
    }
    T rad;
    int i = 0;
    if constexpr (D == 2) {
      Vector<T, 2> center;
      Vector<T, 2> pt;
      for (auto leaf : leafs) {
        center = leaf->getCenter();
        rad = leaf->getRadius();

        pt = center + Vector<T, 2>(-rad, -rad);
        points << pt[0] << " " << pt[1] << " 0 ";
        pt = center + Vector<T, 2>(rad, -rad);
        points << pt[0] << " " << pt[1] << " 0 ";
        pt = center + Vector<T, 2>(rad, rad);
        points << pt[0] << " " << pt[1] << " 0 ";
        pt = center + Vector<T, 2>(-rad, rad);
        points << pt[0] << " " << pt[1] << " 0 \n";

        cells << "4 " << i << " " << i + 1 << " " << i + 2 << " " << i + 3 << " \n";
        i += 4;
      }
    } else if constexpr (D == 3) {
      Vector<T, 3> center;
      Vector<T, 3> pt;
      for (auto leaf : leafs) {
        center = leaf->getCenter();
        rad = leaf->getRadius();

        pt = center + Vector<T, 3>(-rad, -rad, -rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(rad, -rad, -rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(rad, rad, -rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(-rad, rad, -rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(-rad, -rad, rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(rad, -rad, rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(rad, rad, rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " ";
        pt = center + Vector<T, 3>(-rad, rad, rad);
        points << pt[0] << " " << pt[1] << " " << pt[2] << " \n";

        cells << "8 " << i << " " << i + 1 << " " << i + 2 << " " << i + 3 << " " << i + 4
              << " " << i + 5 << " " << i + 6 << " " << i + 7 << " \n";
        i += 8;
      }
    }
    celldataheader << "CELL_DATA " << leafs.size() << std::endl;
    f << points.str() << cells.str() << celldataheader.str();
  }
};

// --------------------old version--------------------
template <typename T, unsigned int D>
class vtkWriterStruPoints {
 private:
  int _Nx;
  int _Ny;
  int _Nz;
  T VoxelSize;
  Vector<T, D> _Min;
  std::string _dirname = "./vtkoutput";
  std::string _filename;
  std::vector<std::vector<int> *> _ScalarIndex;
  std::vector<std::vector<int> *> _VectorIndex;
  std::vector<std::vector<int> *> _FlagIndex;
  std::vector<std::vector<T> *> _Scalars;
  std::vector<std::vector<Vector<T, D>> *> _Vectors;
  std::vector<std::any> _Flags;
  std::vector<std::string> _ScalarNames;
  std::vector<std::string> _VectorNames;
  std::vector<std::string> _FlagNames;

 public:
  vtkWriterStruPoints(std::string filename, T voxelsize, const Vector<T, D> &Min, int Nx,
                      int Ny, int Nz = 1)
      : _filename(filename), VoxelSize(voxelsize), _Min(Min), _Nx(Nx), _Ny(Ny), _Nz(Nz) {}
  void setFilename(std::string filename) { _filename = filename; }
  void setDirname(std::string dirname) { _dirname = dirname; }

  void addScalartoWriteList(std::string name, std::vector<T> *v, std::vector<int> *idx) {
    _ScalarNames.push_back(name);
    _ScalarIndex.push_back(idx);
    _Scalars.push_back(v);
  }
  void addVectortoWriteList(std::string name, std::vector<Vector<T, D>> *v,
                            std::vector<int> *idx) {
    _VectorNames.push_back(name);
    _VectorIndex.push_back(idx);
    _Vectors.push_back(v);
  }
  template <typename U = int>
  void addFlagtoWriteList(std::string name, std::vector<U> *v, std::vector<int> *idx) {
    _FlagNames.push_back(name);
    _FlagIndex.push_back(idx);
    _Flags.push_back(v);
  }
  template <typename U = int>
  void write(int step) {
    DirCreator::Create_Dir(_dirname);
    std::string fullName = "./vtkoutput/" + _filename + std::to_string(step) + ".vtk";
    std::ofstream f(fullName.c_str());
    writeHeader(f);
    // write scalars
    for (int i = 0; i < _Scalars.size(); ++i) {
      writeScalar(f, _ScalarNames[i], *_Scalars[i], *_ScalarIndex[i]);
    }
    // write vectors
    for (int i = 0; i < _Vectors.size(); ++i) {
      writeVector(f, _VectorNames[i], *_Vectors[i], *_VectorIndex[i]);
    }
    // write flags
    for (int i = 0; i < _Flags.size(); ++i) {
      auto flag = std::any_cast<std::vector<U> *>(_Flags[i]);
      WriteFlag(f, _FlagNames[i], *flag, *_FlagIndex[i]);
    }
    f.close();
  }
  template <typename U = int>
  void write() {
    DirCreator::Create_Dir(_dirname);
    std::string fullName = "./vtkoutput/" + _filename + ".vtk";
    std::ofstream f(fullName.c_str());
    writeHeader(f);
    // write scalars
    for (int i = 0; i < _Scalars.size(); ++i) {
      writeScalar(f, _ScalarNames[i], *_Scalars[i], *_ScalarIndex[i]);
    }
    // write vectors
    for (int i = 0; i < _Vectors.size(); ++i) {
      writeVector(f, _VectorNames[i], *_Vectors[i], *_VectorIndex[i]);
    }
    // write flags
    for (int i = 0; i < _Flags.size(); ++i) {
      auto flag = std::any_cast<std::vector<U> *>(_Flags[i]);
      WriteFlag(f, _FlagNames[i], *flag, *_FlagIndex[i]);
    }
    f.close();
  }
  void writeHeader(std::ofstream &f) {
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "Voxels" << std::endl;
    f << "ASCII" << std::endl;

    f << "DATASET STRUCTURED_POINTS" << std::endl;
    f << "DIMENSIONS " << _Nx << " " << _Ny << " " << _Nz << std::endl;
    if constexpr (D == 2) {
      f << "ORIGIN " << _Min[0] + VoxelSize * T(0.5) << " "
        << _Min[1] + VoxelSize * T(0.5) << " 0" << std::endl;
      f << "ASPECT_RATIO " << VoxelSize << " " << VoxelSize << " 1" << std::endl;
    } else {
      f << "ORIGIN " << _Min[0] + VoxelSize * T(0.5) << " "
        << _Min[1] + VoxelSize * T(0.5) << " " << _Min[2] + VoxelSize * T(0.5)
        << std::endl;
      f << "ASPECT_RATIO " << VoxelSize << " " << VoxelSize << " " << VoxelSize
        << std::endl;
    }
    f << "POINT_DATA " << _Nx * _Ny * _Nz << std::endl;
  }

  void writeScalar(std::ofstream &f, std::string name, const std::vector<T> &v,
                   const std::vector<int> &idx) {
    std::stringstream ss;
    if constexpr (std::is_same<T, double>::value) {
      ss << "SCALARS " << name << " double" << std::endl;
    } else if constexpr (std::is_same<T, float>::value) {
      ss << "SCALARS " << name << " float" << std::endl;
    } else if constexpr (std::is_same<T, int>::value) {
      ss << "SCALARS " << name << " int" << std::endl;
    }
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
      if (idx[i] != -1)
        ss << v[idx[i]] << " ";
      else
        ss << 0 << " ";
    }
    ss << std::endl;
    f << ss.str();
  }
  template <typename U = int>
  void WriteFlag(std::ofstream &f, std::string name, const std::vector<U> &v,
                 const std::vector<int> &idx) {
    std::stringstream ss;
    ss << "SCALARS " << name << " int" << std::endl;
    ss << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
      if (idx[i] != -1)
        ss << static_cast<int>(v[idx[i]]) << " ";
      else
        ss << -2 << " ";
    }
    ss << std::endl;
    f << ss.str();
  }

  void writeVector(std::ofstream &f, std::string name, const std::vector<Vector<T, D>> &v,
                   const std::vector<int> &idx) {
    std::stringstream ss;
    if constexpr (std::is_same<T, double>::value) {
      ss << "VECTORS " << name << " double" << std::endl;
    } else if constexpr (std::is_same<T, float>::value) {
      ss << "VECTORS " << name << " float" << std::endl;
    } else if constexpr (std::is_same<T, int>::value) {
      ss << "VECTORS " << name << " int" << std::endl;
    }
    for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
      if (idx[i] != -1) {
        if constexpr (D == 2) {
          ss << v[idx[i]][0] << " " << v[idx[i]][1] << " " << 0 << " ";
        } else {
          ss << v[idx[i]][0] << " " << v[idx[i]][1] << " " << v[idx[i]][2] << " ";
        }
      } else {
        ss << 0 << " " << 0 << " " << 0 << " ";
      }
    }
    ss << std::endl;
    f << ss.str();
  }
  void closeStream(std::ofstream &f) { f.close(); }
};