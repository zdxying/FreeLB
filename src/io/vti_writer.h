/* This file is part of FreeLB
 * 
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 * 
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

// vti_writer.h

#pragma once

#include "data_struct/Vector.h"
#include "data_struct/lattice.h"
#include "io/base64.h"
#include "io/basic_writer.h"

// .pvd - .vtm - .vti
// https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#imagedata
/*
<VTKFile type="ImageData" ...>
  <ImageData WholeExtent="x1 x2 y1 y2 z1 z2"
   Origin="x0 y0 z0" Spacing="dx dy dz">
   <Piece Extent="x1 x2 y1 y2 z1 z2">
      <PointData>...</PointData>
      <CellData>...</CellData>
   </Piece>
   </ImageData>
</VTKFile>
*/

// from The VTK User’s Guide 11th Edition P487
// type of DataArray:
// type — The data type of a single component of the array.
// This is one of
// Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32, Float64.

namespace vtiwriter {

class AbstractWriter {
 public:
  virtual void write(const std::string &fName) = 0;
  virtual void writeBinary(const std::string &fName) = 0;
};

class AbstWriterSet {
 public:
  virtual AbstractWriter *getWriter(int i) = 0;
};

// manager of generating a single vti file
template <typename T, unsigned int D>
class vtiManager {
 private:
  std::string _dirname = "./vtkoutput/vtidata/";
  std::string _filename;
  int _id;

  // image info
  T VoxelSize;
  Vector<T, D> _Min;
  // _Ext = Mesh - Vector<int, D>{1}
  Vector<int, D> _Ext;

#ifdef MPI_ENABLED
  int _rank;
#endif

  std::vector<AbstractWriter *> _writers;

 public:
  vtiManager(std::string filename, int id, T voxsize, Vector<T, D> min, Vector<int, D> ext)
      : _filename(filename), _id(id), VoxelSize(voxsize), _Min(min), _Ext(ext) {
    DirCreator::MPI_Create_Dir(_dirname);
  }
  vtiManager(std::string filename, T voxsize, Vector<T, D> min, Vector<int, D> ext)
      : _filename(filename), _id(-1), VoxelSize(voxsize), _Min(min), _Ext(ext) {
    DirCreator::MPI_Create_Dir(_dirname);
  }
  // extended block should be used
  vtiManager(std::string filename, int id, const BasicBlock<T, D> &block)
      : _filename(filename), _id(id), VoxelSize(block.getVoxelSize()), _Min(block.getMinCenter()),
        _Ext(block.getMesh() - Vector<int, D>{1}) {
    DirCreator::MPI_Create_Dir(_dirname);
  }

  int getBlockId() { return _id; }

  void addWriter(AbstractWriter *writer) { _writers.push_back(writer); }
  template <typename... Args>
  void addWriter(AbstractWriter *writer, Args... args) {
    addWriter(writer);
    addWriter(args...);
  }
  std::string getFileName() {
    std::string fName = _filename + "_B" + std::to_string(_id) + ".vti";
    return fName;
  }
  std::string getFileName(int step) {
    std::string fName =
      _filename + "_T" + std::to_string(step) + "_B" + std::to_string(_id) + ".vti";
    return fName;
  }
  void Write() {
    std::string fullName = _dirname + _filename;
    if (_id == -1) {
      fullName += ".vti";
    } else {
      fullName += "_B" + std::to_string(_id) + ".vti";
    }
    // std::string fullName = _dirname + _filename + "_B" + std::to_string(_id) + ".vti";
    vtiHeader(fullName);
    for (AbstractWriter *writer : _writers) {
      writer->write(fullName);
    }
    vtiEnd(fullName);
  }
  void Write(int step) {
    std::string fullName = _dirname + _filename;
    if (_id == -1) {
      fullName += "_T" + std::to_string(step) + ".vti";
    } else {
      fullName += "_T" + std::to_string(step) + "_B" + std::to_string(_id) + ".vti";
    }
    // std::string fullName =
    // _dirname + _filename + "_T" + std::to_string(step) + "_B" + std::to_string(_id) + ".vti";
    vtiHeader(fullName);
    for (AbstractWriter *writer : _writers) {
      writer->write(fullName);
    }
    vtiEnd(fullName);
  }
  void WriteBinary() {
    std::string fullName = _dirname + _filename;
    if (_id == -1) {
      fullName += ".vti";
    } else {
      fullName += "_B" + std::to_string(_id) + ".vti";
    }
    vtiHeader(fullName);

    for (AbstractWriter *writer : _writers) {
      writer->writeBinary(fullName);
    }
    vtiEnd(fullName);
  }
  void WriteBinary(int step) {
    std::string fullName = _dirname + _filename;
    if (_id == -1) {
      fullName += "_T" + std::to_string(step) + ".vti";
    } else {
      fullName += "_T" + std::to_string(step) + "_B" + std::to_string(_id) + ".vti";
    }
    // std::string fullName =
    // _dirname + _filename + "_T" + std::to_string(step) + "_B" + std::to_string(_id) + ".vti";
    vtiHeader(fullName);
    for (AbstractWriter *writer : _writers) {
      writer->writeBinary(fullName);
    }
    vtiEnd(fullName);
  }
  void vtiHeader(const std::string &fName) {
    std::ofstream f(fName);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" "
      << "byte_order=\"LittleEndian\">\n";
    if constexpr (D == 2) {
      f << "<ImageData WholeExtent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " "
        << 0 << " " << 0 << "\" Origin=\"" << _Min[0] << " " << _Min[1] << " " << 0
        << "\" Spacing=\"" << VoxelSize << " " << VoxelSize << " " << VoxelSize << "\">\n";
      f << "<Piece Extent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " " << 0
        << " " << 0 << "\">\n";
    } else if constexpr (D == 3) {
      f << "<ImageData WholeExtent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " "
        << 0 << " " << _Ext[2] << "\" Origin=\"" << _Min[0] << " " << _Min[1] << " " << _Min[2]
        << "\" Spacing=\"" << VoxelSize << " " << VoxelSize << " " << VoxelSize << "\">\n";
      f << "<Piece Extent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " " << 0
        << " " << _Ext[2] << "\">\n";
    }
    f << "<PointData>\n";
    f.close();
  }
  void vtiEnd(const std::string &fName) {
    std::ofstream f(fName, std::ios::app);
    f << "</PointData>\n";
    f << "</Piece>\n";
    f << "</ImageData>\n";
    f << "</VTKFile>\n";
    f.close();
  }
};

template <template <typename> class ArrayType, typename T>
class ScalerWriter : public AbstractWriter {
 private:
  std::string varname;
  // field data
  const ArrayType<T> &Array;
  std::size_t Size;

 public:
  ScalerWriter(std::string name, const ArrayType<T> &f) : varname(name), Array(f), Size(f.size()) {}
  void write(const std::string &fName) override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    // getVTKTypeString<T>(type);
    f << "<DataArray type=\""
      << "Float32"
      << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << 1 << "\">\n";
    for (int i = 0; i < Size; ++i) {
      f << static_cast<float>(Array[i]) << " ";
    }
    f << "\n</DataArray>\n";
    f.close();
  }
  void writeBinary(const std::string &fName) override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<T>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" "
      << "NumberOfComponents=\"" << 1 << "\">\n";
    f.close();

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(Size * sizeof(T));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<T> Encoder(fb, Size);
    Encoder.encode(Array.getdata(), Size);
    fb.close();

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }
};

template <typename T, unsigned int D>
class VectorWriter : public AbstractWriter {
 private:
  std::string varname;
  // field data
  const GenericArray<Vector<T, D>> &Array;
  std::size_t Size;

 public:
  VectorWriter(std::string name, const GenericArray<Vector<T, D>> &f)
      : varname(name), Array(f), Size(f.size()) {}
  void write(const std::string &fName) override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    // getVTKTypeString<T>(type);
    f << "<DataArray type=\""
      << "Float32"
      << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << D << "\">\n";
    for (int i = 0; i < Size; ++i) {
      if constexpr (D == 2) {
        f << static_cast<float>(Array[i][0]) << " " << static_cast<float>(Array[i][1]) << " ";
      } else if constexpr (D == 3) {
        f << static_cast<float>(Array[i][0]) << " " << static_cast<float>(Array[i][1]) << " "
          << static_cast<float>(Array[i][2]) << " ";
      }
    }
    f << "\n</DataArray>\n";
    f.close();
  }
  void writeBinary(const std::string &fName) override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<T>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" "
      << "NumberOfComponents=\"" << D << "\">\n";
    f.close();

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t DSize = std::size_t(D * Size);
    std::size_t binarySize = std::size_t(DSize * sizeof(T));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);

    // get data array for encoding
    T *data = new T[DSize];
    if constexpr (D == 2 || D == 3) {
      int j = 0;
      for (int i = 0; i < Size; ++i) {
        data[j] = Array[i][0];
        data[j + 1] = Array[i][1];
        if constexpr (D == 3) {
          data[j + 2] = Array[i][2];
        }
        j += D;
      }
    } else {
      int j = 0;
      for (int i = 0; i < Size; i += D) {
        for (unsigned int k = 0; k < D; ++k) {
          data[j + k] = Array[i][k];
        }
        j += D;
      }
    }
    // writes the data
    Base64Encoder<T> Encoder(fb, DSize);
    Encoder.encode(data, DSize);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }
};

template <template <typename> class ArrayType, typename T, unsigned int D>
class VectorSOAWriter : public AbstractWriter {
 private:
  std::string varname;
  // field data
  const GenericField<ArrayType<T>, D> &Field;
  std::size_t Size;

 public:
  VectorSOAWriter(std::string name, const GenericField<ArrayType<T>, D> &f)
      : varname(name), Field(f), Size(f.getField(0).size()) {}
  void write(const std::string &fName) override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    // getVTKTypeString<T>(type);
    f << "<DataArray type=\""
      << "Float32"
      << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << D << "\">\n";
    for (int i = 0; i < Size; ++i) {
      for (unsigned int k = 0; k < D; ++k) {
        f << static_cast<float>(Field.template get<k>(i)) << " ";
      }
    }
    f << "\n</DataArray>\n";
    f.close();
  }
  void writeBinary(const std::string &fName) override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<T>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" "
      << "NumberOfComponents=\"" << D << "\">\n";
    f.close();

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t DSize = std::size_t(D * Size);
    std::size_t binarySize = std::size_t(DSize * sizeof(T));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);

    // get data array for encoding
    T *data = new T[DSize];
    int j = 0;
    for (int i = 0; i < Size; i += D) {
      for (unsigned int k = 0; k < D; ++k) {
        data[j + k] = Field.get(i, k);
      }
      j += D;
    }
    // writes the data
    Base64Encoder<T> Encoder(fb, DSize);
    Encoder.encode(data, DSize);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }
};


}  // namespace vtiwriter

// vti writer with overlap correction

// non-overlapped vti writer
namespace vtino {

class AbstractWriter {
 public:
  // virtual void write(const std::string& fName) = 0;
  virtual void writeBinary(const std::string &fName, int Overlap) = 0;
};

class AbstWriterSet {
 public:
  virtual AbstractWriter *getWriter(int i) = 0;
};

// manager of generating a single vti file
template <typename T, unsigned int D>
class vtiManager {
 private:
  std::string _dirname = "./vtkoutput/vtidata/";
  std::string _filename;
  int _id;

  // image info
  T VoxelSize;
  Vector<T, D> _Min;
  // _Ext = Mesh - Vector<int, D>{1}
  Vector<int, D> _Ext;
  // overlap
  int _Overlap;

#ifdef MPI_ENABLED
  int _rank;
#endif

  std::vector<AbstractWriter *> _writers;

 public:
  vtiManager(std::string filename, int id, T voxsize, Vector<T, D> min, Vector<int, D> ext,
             int Overlap)
      : _filename(filename), _id(id), VoxelSize(voxsize), _Overlap(Overlap),
        _Min(min + Vector<T, D>{Overlap * voxsize}), _Ext(ext - Vector<int, D>{2 * Overlap}) {
    DirCreator::MPI_Create_Dir(_dirname);
  }
  // extended block should be used
  vtiManager(std::string filename, int id, const BasicBlock<T, D> &block, int Overlap)
      : _filename(filename), _id(id), VoxelSize(block.getVoxelSize()), _Overlap(Overlap),
        _Min(block.getMinCenter() + Vector<T, D>{Overlap * block.getVoxelSize()}),
        _Ext(block.getMesh() - Vector<int, D>{1} - Vector<int, D>{2 * Overlap}) {
    DirCreator::MPI_Create_Dir(_dirname);
  }

  int getBlockId() { return _id; }

  void addWriter(AbstractWriter *writer) { _writers.push_back(writer); }
  template <typename... Args>
  void addWriter(AbstractWriter *writer, Args... args) {
    addWriter(writer);
    addWriter(args...);
  }
  std::string getFileName() {
    std::string fName = _filename + "_B" + std::to_string(_id) + ".vti";
    return fName;
  }
  std::string getFileName(int step) {
    std::string fName =
      _filename + "_T" + std::to_string(step) + "_B" + std::to_string(_id) + ".vti";
    return fName;
  }
  void WriteBinary() {
    std::string fullName = _dirname + _filename;
    if (_id == -1) {
      fullName += ".vti";
    } else {
      fullName += "_B" + std::to_string(_id) + ".vti";
    }
    vtiHeader(fullName);
    for (AbstractWriter *writer : _writers) {
      writer->writeBinary(fullName, _Overlap);
    }
    vtiEnd(fullName);
  }
  void WriteBinary(int step) {
    std::string fullName = _dirname + _filename;
    if (_id == -1) {
      fullName += "_T" + std::to_string(step) + ".vti";
    } else {
      fullName += "_T" + std::to_string(step) + "_B" + std::to_string(_id) + ".vti";
    }
    vtiHeader(fullName);
    for (AbstractWriter *writer : _writers) {
      writer->writeBinary(fullName, _Overlap);
    }
    vtiEnd(fullName);
  }
  void vtiHeader(const std::string &fName) {
    std::ofstream f(fName);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" "
      << "byte_order=\"LittleEndian\">\n";
    if constexpr (D == 2) {
      f << "<ImageData WholeExtent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " "
        << 0 << " " << 0 << "\" Origin=\"" << _Min[0] << " " << _Min[1] << " " << 0
        << "\" Spacing=\"" << VoxelSize << " " << VoxelSize << " " << VoxelSize << "\">\n";
      f << "<Piece Extent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " " << 0
        << " " << 0 << "\">\n";
    } else if constexpr (D == 3) {
      f << "<ImageData WholeExtent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " "
        << 0 << " " << _Ext[2] << "\" Origin=\"" << _Min[0] << " " << _Min[1] << " " << _Min[2]
        << "\" Spacing=\"" << VoxelSize << " " << VoxelSize << " " << VoxelSize << "\">\n";
      f << "<Piece Extent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " " << 0
        << " " << _Ext[2] << "\">\n";
    }
    f << "<PointData>\n";
    f.close();
  }
  void vtiEnd(const std::string &fName) {
    std::ofstream f(fName, std::ios::app);
    f << "</PointData>\n";
    f << "</Piece>\n";
    f << "</ImageData>\n";
    f << "</VTKFile>\n";
    f.close();
  }
};

template <template <typename> class ArrayType, typename T, unsigned int D>
class ScalerWriter : public AbstractWriter {
 private:
  std::string varname;
  // field data
  const ArrayType<T> &Array;
  // mesh info
  Vector<int, D> Mesh;

 public:
  ScalerWriter(std::string name, const ArrayType<T> &f, Vector<int, D> mesh)
      : varname(name), Array(f), Mesh(mesh) {}

  void writeBinary(const std::string &fName, int Overlap) override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<T>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" "
      << "NumberOfComponents=\"" << 1 << "\">\n";
    f.close();

    // get submesh size
    std::size_t Size;
    Vector<int, D> SubMesh = Mesh - Vector<int, D>{2 * Overlap};
    if constexpr (D == 2) {
      Size = SubMesh[0] * SubMesh[1];
    } else if constexpr (D == 3) {
      Size = SubMesh[0] * SubMesh[1] * SubMesh[2];
    }
    // get data array for encoding
    T *data = new T[Size];
    std::size_t id = 0;
    if constexpr (D == 2) {
      for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
        // copy data from {Array[Overlap][j] to Array[Mesh[0]-Overlap-1][j]} to data
        std::size_t start = Overlap + j * Mesh[0];
        std::copy(Array.getdataPtr(start), Array.getdataPtr(start + SubMesh[0]), data + id);
        id += SubMesh[0];
      }
    } else if constexpr (D == 3) {
      int XY = Mesh[0] * Mesh[1];
      for (int k = Overlap; k < Mesh[2] - Overlap; ++k) {
        for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
          // copy data from {Array[Overlap][j][k] to Array[Mesh[0]-Overlap-1][j][k]} to data
          int start = Overlap + j * Mesh[0] + k * XY;
          std::copy(Array.getdataPtr(start), Array.getdataPtr(start + SubMesh[0]), data + id);
          id += SubMesh[0];
        }
      }
    }

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(Size * sizeof(T));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<T> Encoder(fb, Size);
    Encoder.encode(data, Size);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }
};

template <typename T, unsigned int D, unsigned int arrD>
class VectorWriter : public AbstractWriter {
 private:
  std::string varname;
  // field data
  const GenericArray<Vector<T, arrD>> &Array;
  // mesh info
  Vector<int, D> Mesh;

 public:
  VectorWriter(std::string name, const GenericArray<Vector<T, arrD>> &f, Vector<int, D> mesh)
      : varname(name), Array(f), Mesh(mesh) {}

  void writeBinary(const std::string &fName, int Overlap) override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<T>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" "
      << "NumberOfComponents=\"" << arrD << "\">\n";
    f.close();

    // get submesh size
    std::size_t Size;
    Vector<int, D> SubMesh = Mesh - Vector<int, D>{2 * Overlap};
    if constexpr (D == 2) {
      Size = SubMesh[0] * SubMesh[1];
    } else if constexpr (D == 3) {
      Size = SubMesh[0] * SubMesh[1] * SubMesh[2];
    }
    std::size_t arrDSize = std::size_t(arrD * Size);

    // get data array for encoding
    T *data = new T[arrDSize];
    std::size_t id = 0;
    if constexpr (D == 2) {
      for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
        for (int i = Overlap; i < Mesh[0] - Overlap; ++i) {
          std::size_t arrId = i + j * Mesh[0];
          std::copy(Array[arrId].data(), Array[arrId].data() + arrD, data + id);
          id += arrD;
        }
      }
    } else if constexpr (D == 3) {
      int XY = Mesh[0] * Mesh[1];
      for (int k = Overlap; k < Mesh[2] - Overlap; ++k) {
        for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
          for (int i = Overlap; i < Mesh[0] - Overlap; ++i) {
            std::size_t arrId = i + j * Mesh[0] + k * XY;
            std::copy(Array[arrId].data(), Array[arrId].data() + arrD, data + id);
            id += arrD;
          }
        }
      }
    }

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(arrDSize * sizeof(T));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<T> Encoder(fb, arrDSize);
    Encoder.encode(data, arrDSize);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }
};

template <template <typename> class ArrayType, typename T, unsigned int D, unsigned int arrD>
class VectorSOAWriter : public AbstractWriter {
 private:
  std::string varname;
  // field data
  const GenericField<ArrayType<T>, arrD> &Field;
  // mesh info
  Vector<int, D> Mesh;

 public:
  VectorSOAWriter(std::string name, const GenericField<ArrayType<T>, arrD> &f, Vector<int, D> mesh)
      : varname(name), Field(f), Mesh(mesh) {}

  void writeBinary(const std::string &fName, int Overlap) override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<T>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" "
      << "NumberOfComponents=\"" << arrD << "\">\n";
    f.close();

    // get submesh size
    std::size_t Size;
    Vector<int, D> SubMesh = Mesh - Vector<int, D>{2 * Overlap};
    if constexpr (D == 2) {
      Size = SubMesh[0] * SubMesh[1];
    } else if constexpr (D == 3) {
      Size = SubMesh[0] * SubMesh[1] * SubMesh[2];
    }
    std::size_t arrDSize = std::size_t(arrD * Size);

    // get data array for encoding
    T *data = new T[arrDSize];
    std::size_t id = 0;
    if constexpr (D == 2) {
      for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
        for (int i = Overlap; i < Mesh[0] - Overlap; ++i) {
          std::size_t arrId = i + j * Mesh[0];
          for (unsigned int vecid = 0; vecid < arrD; ++vecid) {
            data[id + vecid] = Field.get(arrId, vecid);
          }
          id += arrD;
        }
      }
    } else if constexpr (D == 3) {
      int XY = Mesh[0] * Mesh[1];
      for (int k = Overlap; k < Mesh[2] - Overlap; ++k) {
        for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
          for (int i = Overlap; i < Mesh[0] - Overlap; ++i) {
            std::size_t arrId = i + j * Mesh[0] + k * XY;
            for (unsigned int vecid = 0; vecid < arrD; ++vecid) {
              data[id + vecid] = Field.get(arrId, vecid);
            }
            id += arrD;
          }
        }
      }
    }

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(arrDSize * sizeof(T));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);

    Base64Encoder<T> Encoder(fb, arrDSize);
    Encoder.encode(data, arrDSize);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }
};

}  // namespace vtino