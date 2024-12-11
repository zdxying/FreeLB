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

// vti_writer.h

#pragma once

// store any callable object and used in unit converter when writing vti files
#include <functional>

#include "data_struct/Vector.h"
#include "data_struct/block_lattice.h"
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
  virtual void write(const std::string &fName) const = 0;
  virtual void writeBinary(const std::string &fName) const = 0;
};

class AbstWriterSet {
 public:
  virtual const AbstractWriter &getWriter(int i) const = 0;
};

// manager of generating a single vti file
template <typename T, unsigned int Dim>
class vtiManager {
 private:
  std::string _dirname = "./vtkoutput/vtidata/";
  std::string _filename;
  int _id;

  // image info
  T VoxelSize;
  Vector<T, Dim> _Min;
  // _Ext = Mesh - Vector<int, Dim>{1}
  Vector<int, Dim> _Ext;

#ifdef MPI_ENABLED
  int _rank;
#endif

  std::vector<const AbstractWriter *> _writers;

 public:
  vtiManager(std::string filename, int id, T voxsize, Vector<T, Dim> min,
             Vector<int, Dim> ext)
      : _filename(filename), _id(id), VoxelSize(voxsize), _Min(min), _Ext(ext) {
    static_assert(Dim == 2 || Dim == 3, "Error: Dimension is not supported!");
    DirCreator::Create_Dir(_dirname);
  }
  vtiManager(std::string filename, T voxsize, Vector<T, Dim> min, Vector<int, Dim> ext)
      : _filename(filename), _id(-1), VoxelSize(voxsize), _Min(min), _Ext(ext) {
    static_assert(Dim == 2 || Dim == 3, "Error: Dimension is not supported!");
    DirCreator::Create_Dir(_dirname);
  }
  // extended block should be used
  vtiManager(std::string filename, int id, const BasicBlock<T, Dim> &block)
      : _filename(filename), _id(id), VoxelSize(block.getVoxelSize()),
        _Min(block.getMinCenter()), _Ext(block.getMesh() - Vector<int, Dim>{1}) {
    static_assert(Dim == 2 || Dim == 3, "Error: Dimension is not supported!");
    DirCreator::Create_Dir(_dirname);
  }

  int getvtiBlockId() { return _id; }

  void addWriter(const AbstractWriter &writer) { _writers.push_back(&writer); }
  template <typename... Args>
  void addWriter(const AbstractWriter &writer, const Args &...args) {
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
    for (const AbstractWriter *writer : _writers) {
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
    // _dirname + _filename + "_T" + std::to_string(step) + "_B" + std::to_string(_id) +
    // ".vti";
    vtiHeader(fullName);
    for (const AbstractWriter *writer : _writers) {
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

    for (const AbstractWriter *writer : _writers) {
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
    // _dirname + _filename + "_T" + std::to_string(step) + "_B" + std::to_string(_id) +
    // ".vti";
    vtiHeader(fullName);
    for (const AbstractWriter *writer : _writers) {
      writer->writeBinary(fullName);
    }
    vtiEnd(fullName);
  }
  void vtiHeader(const std::string &fName) {
    std::ofstream f(fName);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" "
      << "byte_order=\"LittleEndian\">\n";
    if constexpr (Dim == 2) {
      f << "<ImageData WholeExtent=\"" << 0 << " " << _Ext[0] << " " << 0 << " "
        << _Ext[1] << " " << 0 << " " << 0 << "\" Origin=\"" << _Min[0] << " " << _Min[1]
        << " " << 0 << "\" Spacing=\"" << VoxelSize << " " << VoxelSize << " "
        << VoxelSize << "\">\n";
      f << "<Piece Extent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " "
        << 0 << " " << 0 << "\">\n";
    } else if constexpr (Dim == 3) {
      f << "<ImageData WholeExtent=\"" << 0 << " " << _Ext[0] << " " << 0 << " "
        << _Ext[1] << " " << 0 << " " << _Ext[2] << "\" Origin=\"" << _Min[0] << " "
        << _Min[1] << " " << _Min[2] << "\" Spacing=\"" << VoxelSize << " " << VoxelSize
        << " " << VoxelSize << "\">\n";
      f << "<Piece Extent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " "
        << 0 << " " << _Ext[2] << "\">\n";
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

template <typename ArrayType>
class ScalarWriter : public AbstractWriter {
 private:
  std::string varname;
  // field data
  const ArrayType &Array;
  std::size_t Size;

 public:
  using datatype = typename ArrayType::value_type;

  ScalarWriter(std::string name, const ArrayType &f)
      : varname(name), Array(f), Size(f.size()) {}
  void write(const std::string &fName) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    // getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << "Float32" << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << 1 << "\">\n";
    for (std::size_t i = 0; i < Size; ++i) {
      f << static_cast<float>(Array[i]) << " ";
    }
    f << "\n</DataArray>\n";
    f.close();
  }
  void writeBinary(const std::string &fName) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" " << "NumberOfComponents=\"" << 1
      << "\">\n";
    f.close();

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(Size * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<datatype> Encoder(fb, Size);
    Encoder.encode(Array.getdataPtr(), Size);
    fb.close();

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }
};

template <typename ArrayType>
class VectorWriter : public AbstractWriter {
 public:
  using vectortype = typename ArrayType::value_type;
  using datatype = typename vectortype::value_type;
  static constexpr unsigned int D = vectortype::vector_dim;

 private:
  std::string varname;
  // field data
  const ArrayType &Array;
  std::size_t Size;

 public:
  VectorWriter(std::string name, const ArrayType &f)
      : varname(name), Array(f), Size(f.size()) {}
  void write(const std::string &fName) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    // getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << "Float32" << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << D << "\">\n";
    for (std::size_t i = 0; i < Size; ++i) {
      if constexpr (D == 2) {
        f << static_cast<float>(Array[i][0]) << " " << static_cast<float>(Array[i][1])
          << " ";
      } else if constexpr (D == 3) {
        f << static_cast<float>(Array[i][0]) << " " << static_cast<float>(Array[i][1])
          << " " << static_cast<float>(Array[i][2]) << " ";
      }
    }
    f << "\n</DataArray>\n";
    f.close();
  }
  void writeBinary(const std::string &fName) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" " << "NumberOfComponents=\"" << D
      << "\">\n";
    f.close();

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t DSize = std::size_t(D * Size);
    std::size_t binarySize = std::size_t(DSize * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);

    // get data array for encoding
    datatype *data = new datatype[DSize];
    if constexpr (D == 2 || D == 3) {
      std::size_t j = 0;
      for (std::size_t i = 0; i < Size; ++i) {
        data[j] = Array[i][0];
        data[j + 1] = Array[i][1];
        if constexpr (D == 3) {
          data[j + 2] = Array[i][2];
        }
        j += D;
      }
    } else {
      std::size_t j = 0;
      for (std::size_t i = 0; i < Size; i += D) {
        for (unsigned int k = 0; k < D; ++k) {
          data[j + k] = Array[i][k];
        }
        j += D;
      }
    }
    // writes the data
    Base64Encoder<datatype> Encoder(fb, DSize);
    Encoder.encode(data, DSize);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }
};

template <typename ArrayType, unsigned int D>
class VectorSOAWriter : public AbstractWriter {
 private:
  std::string varname;
  // field data
  const GenericArrayField<ArrayType, D> &Field;
  std::size_t Size;

 public:
  using datatype = typename ArrayType::value_type;

  VectorSOAWriter(std::string name, const GenericArrayField<ArrayType, D> &f)
      : varname(name), Field(f), Size(f.getField(0).size()) {}
  void write(const std::string &fName) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    // getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << "Float32" << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << D << "\">\n";
    for (std::size_t i = 0; i < Size; ++i) {
      for (unsigned int k = 0; k < D; ++k) {
        f << static_cast<float>(Field.template get<k>(i)) << " ";
      }
    }
    f << "\n</DataArray>\n";
    f.close();
  }
  void writeBinary(const std::string &fName) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" " << "NumberOfComponents=\"" << D
      << "\">\n";
    f.close();

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t DSize = std::size_t(D * Size);
    std::size_t binarySize = std::size_t(DSize * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);

    // get data array for encoding
    datatype *data = new datatype[DSize];
    std::size_t j = 0;
    for (std::size_t i = 0; i < Size; i += D) {
      for (unsigned int k = 0; k < D; ++k) {
        data[j + k] = Field.get(i, k);
      }
      j += D;
    }
    // writes the data
    Base64Encoder<datatype> Encoder(fb, DSize);
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
  virtual void writeBinary(const std::string &fName, int Overlap) const = 0;
  virtual void writeConvertedBinary(const std::string &fName, int Overlap) const = 0;
};

class AbstWriterSet {
 public:
  virtual const AbstractWriter &getWriter(int i) const = 0;
};

// manager of generating a single vti file
template <typename T, unsigned int Dim>
class vtiManager {
 private:
  std::string _dirname = "./vtkoutput/vtidata/";
  std::string _filename;
  int _id;

  // image info
  T VoxelSize;
  // overlap
  int _Overlap;
  Vector<T, Dim> _Min;
  // _Ext = Mesh - Vector<int, Dim>{1}
  Vector<int, Dim> _Ext;

#ifdef MPI_ENABLED
  int _rank;
#endif

  std::vector<const AbstractWriter *> _writers;

 public:
  vtiManager(std::string filename, int id, T voxsize, Vector<T, Dim> min,
             Vector<int, Dim> ext, int Overlap)
      : _filename(filename), _id(id), VoxelSize(voxsize), _Overlap(Overlap),
        _Min(min + Vector<T, Dim>{Overlap * voxsize}),
        _Ext(ext - Vector<int, Dim>{2 * Overlap}) {
    static_assert(Dim == 2 || Dim == 3, "Error: Dimension is not supported!");
    DirCreator::Create_Dir(_dirname);
  }
  // extended block should be used
  vtiManager(std::string filename, int id, const BasicBlock<T, Dim> &block, int Overlap)
      : _filename(filename), _id(id), VoxelSize(block.getVoxelSize()), _Overlap(Overlap),
        _Min(block.getMinCenter() + Vector<T, Dim>{Overlap * block.getVoxelSize()}),
        _Ext(block.getMesh() - Vector<int, Dim>{1} - Vector<int, Dim>{2 * Overlap}) {
    static_assert(Dim == 2 || Dim == 3, "Error: Dimension is not supported!");
    DirCreator::Create_Dir(_dirname);
  }

  int getvtiBlockId() { return _id; }

  void addWriter(const AbstractWriter &writer) { _writers.push_back(&writer); }
  template <typename... Args>
  void addWriter(const AbstractWriter &writer, const Args &...args) {
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
    for (const AbstractWriter *writer : _writers) {
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
    for (const AbstractWriter *writer : _writers) {
      writer->writeBinary(fullName, _Overlap);
    }
    vtiEnd(fullName);
  }
  void WriteConvertedBinary() {
    std::string fullName = _dirname + _filename;
    if (_id == -1) {
      fullName += ".vti";
    } else {
      fullName += "_B" + std::to_string(_id) + ".vti";
    }
    vtiHeader(fullName);
    for (const AbstractWriter *writer : _writers) {
      writer->writeConvertedBinary(fullName, _Overlap);
    }
    vtiEnd(fullName);
  }
  void WriteConvertedBinary(int step) {
    std::string fullName = _dirname + _filename;
    if (_id == -1) {
      fullName += "_T" + std::to_string(step) + ".vti";
    } else {
      fullName += "_T" + std::to_string(step) + "_B" + std::to_string(_id) + ".vti";
    }
    vtiHeader(fullName);
    for (const AbstractWriter *writer : _writers) {
      writer->writeConvertedBinary(fullName, _Overlap);
    }
    vtiEnd(fullName);
  }

  void vtiHeader(const std::string &fName) {
    std::ofstream f(fName);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" "
      << "byte_order=\"LittleEndian\">\n";
    if constexpr (Dim == 2) {
      f << "<ImageData WholeExtent=\"" << 0 << " " << _Ext[0] << " " << 0 << " "
        << _Ext[1] << " " << 0 << " " << 0 << "\" Origin=\"" << _Min[0] << " " << _Min[1]
        << " " << 0 << "\" Spacing=\"" << VoxelSize << " " << VoxelSize << " "
        << VoxelSize << "\">\n";
      f << "<Piece Extent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " "
        << 0 << " " << 0 << "\">\n";
    } else if constexpr (Dim == 3) {
      f << "<ImageData WholeExtent=\"" << 0 << " " << _Ext[0] << " " << 0 << " "
        << _Ext[1] << " " << 0 << " " << _Ext[2] << "\" Origin=\"" << _Min[0] << " "
        << _Min[1] << " " << _Min[2] << "\" Spacing=\"" << VoxelSize << " " << VoxelSize
        << " " << VoxelSize << "\">\n";
      f << "<Piece Extent=\"" << 0 << " " << _Ext[0] << " " << 0 << " " << _Ext[1] << " "
        << 0 << " " << _Ext[2] << "\">\n";
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

template <typename ArrayType, unsigned int Dim>
class ScalarWriter : public AbstractWriter {
 public:
  using datatype = typename ArrayType::value_type;

  ScalarWriter(std::string name, const ArrayType &f, Vector<int, Dim> mesh)
      : varname(name), Array(f), Mesh(mesh) {
    static_assert(Dim == 2 || Dim == 3, "Error: Dimension is not supported!");
  }
  // std::bind(&uintConvclass::func, &unitConv, std::placeholders::_1); or 
  // [&unitConv](T x) { return unitConv.func(x); };
  ScalarWriter(std::string name, const ArrayType &f, Vector<int, Dim> mesh, 
  std::function<datatype(datatype)> func)
      : varname(name), Array(f), Mesh(mesh), unitConvert(func) {
    static_assert(Dim == 2 || Dim == 3, "Error: Dimension is not supported!");
  }

  void writeBinary(const std::string &fName, int Overlap) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" " << "NumberOfComponents=\"" << 1
      << "\">\n";
    f.close();

    datatype *data = nullptr;
    std::size_t Size;
    if (Overlap == 0) {
      Size = Array.size();
    } else {
      Vector<int, Dim> SubMesh = Mesh - Vector<int, Dim>{2 * Overlap};
      if constexpr (Dim == 2) {
        Size = SubMesh[0] * SubMesh[1];
      } else if constexpr (Dim == 3) {
        Size = SubMesh[0] * SubMesh[1] * SubMesh[2];
      }
      // data array for encoding
      data = new datatype[Size];
      // copy data
      util::CopyFromFieldArray(Mesh, Overlap, Array, data);
    }

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(Size * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<datatype> Encoder(fb, Size);
    if (Overlap == 0) {
      Encoder.encode(Array.getdataPtr(), Size);
    } else {
      Encoder.encode(data, Size);
      delete[] data;
    }
    fb.close();

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

  void writeConvertedBinary(const std::string &fName, int Overlap) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" " << "NumberOfComponents=\"" << 1
      << "\">\n";
    f.close();

    datatype *data = nullptr;
    std::size_t Size;
    if (Overlap == 0) {
      Size = Array.size();
    } else {
      Vector<int, Dim> SubMesh = Mesh - Vector<int, Dim>{2 * Overlap};
      if constexpr (Dim == 2) {
        Size = SubMesh[0] * SubMesh[1];
      } else if constexpr (Dim == 3) {
        Size = SubMesh[0] * SubMesh[1] * SubMesh[2];
      }
      // data array for encoding
      data = new datatype[Size];
      // copy data
      util::CopyFromFieldArray(Mesh, Overlap, Array, data, unitConvert);
    }

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(Size * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<datatype> Encoder(fb, Size);
    if (Overlap == 0) {
      Encoder.encode(Array.getdataPtr(), Size);
    } else {
      Encoder.encode(data, Size);
      delete[] data;
    }
    fb.close();

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

 private:
  std::string varname;
  // field data
  const ArrayType &Array;
  // mesh info
  Vector<int, Dim> Mesh;
  // unit convert function pointer
  std::function<datatype(datatype)> unitConvert;
};

template <typename ArrayType, unsigned int Dim>
class VectorWriter : public AbstractWriter {
  public:
  using vectortype = typename ArrayType::value_type;
  using datatype = typename vectortype::value_type;
  static constexpr unsigned int D = vectortype::vector_dim;

  VectorWriter(std::string name, const ArrayType &f, Vector<int, Dim> mesh)
      : varname(name), Array(f), Mesh(mesh) {
    static_assert(Dim == 2 || Dim == 3, "Error: Dimension is not supported!");
  }
  // std::bind(&uintConvclass::func, &unitConv, std::placeholders::_1); or 
  // [&unitConv](T x) { return unitConv.func(x); };
  VectorWriter(std::string name, const ArrayType &f, Vector<int, Dim> mesh, 
  std::function<vectortype(vectortype)> func)
      : varname(name), Array(f), Mesh(mesh), unitConvert(func) {
    static_assert(Dim == 2 || Dim == 3, "Error: Dimension is not supported!");
  }

  void writeBinary(const std::string &fName, int Overlap) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" " << "NumberOfComponents=\"" << D
      << "\">\n";
    f.close();

    Vector<datatype, D> *data = nullptr;
    std::size_t Size;
    std::size_t arrDSize;
    if (Overlap == 0) {
      Size = Array.size();
      arrDSize = D * Size;
    } else {
      Vector<int, Dim> SubMesh = Mesh - Vector<int, Dim>{2 * Overlap};
      if constexpr (Dim == 2) {
        Size = SubMesh[0] * SubMesh[1];
      } else if constexpr (Dim == 3) {
        Size = SubMesh[0] * SubMesh[1] * SubMesh[2];
      }
      arrDSize = D * Size;
      // data array for encoding
      data = new Vector<datatype, D>[arrDSize];
      // copy data
      util::CopyFromFieldArray(Mesh, Overlap, Array, data);
    }

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(arrDSize * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<datatype> Encoder(fb, arrDSize);
    if (Overlap == 0) {
      Encoder.encode(Array.getdataPtr()->data(), arrDSize);
    } else {
      Encoder.encode(data->data(), arrDSize);
      delete[] data;
    }
    fb.close();

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

  void writeConvertedBinary(const std::string &fName, int Overlap) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" " << "NumberOfComponents=\"" << D
      << "\">\n";
    f.close();

    Vector<datatype, D> *data = nullptr;
    std::size_t Size;
    std::size_t arrDSize;
    if (Overlap == 0) {
      Size = Array.size();
      arrDSize = D * Size;
    } else {
      Vector<int, Dim> SubMesh = Mesh - Vector<int, Dim>{2 * Overlap};
      if constexpr (Dim == 2) {
        Size = SubMesh[0] * SubMesh[1];
      } else if constexpr (Dim == 3) {
        Size = SubMesh[0] * SubMesh[1] * SubMesh[2];
      }
      arrDSize = D * Size;
      // data array for encoding
      data = new Vector<datatype, D>[arrDSize];
      // copy data
      util::CopyFromFieldArray(Mesh, Overlap, Array, data, unitConvert);
    }

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(arrDSize * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<datatype> Encoder(fb, arrDSize);
    if (Overlap == 0) {
      Encoder.encode(Array.getdataPtr()->data(), arrDSize);
    } else {
      Encoder.encode(data->data(), arrDSize);
      delete[] data;
    }
    fb.close();

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

 private:
  std::string varname;
  // field data
  const ArrayType &Array;
  // mesh info
  Vector<int, Dim> Mesh;
  // unit convert function pointer
  std::function<vectortype(vectortype)> unitConvert;
};

// XML does not support the SOA format
// we may change to HDF5 format when writing SOA data
// TODO: implement the SOA writer
template <typename ArrayType, unsigned int Dim, unsigned int D>
class VectorSOAWriter : public AbstractWriter {
 public:
  using datatype = typename ArrayType::value_type;

  VectorSOAWriter(std::string name, const GenericArrayField<ArrayType, D> &f, Vector<int, Dim> mesh)
      : varname(name), Field(f), Mesh(mesh) {}
  // std::bind(&uintConvclass::func, &unitConv, std::placeholders::_1); or 
  // [&unitConv](T x) { return unitConv.func(x); };
  VectorSOAWriter(std::string name, const GenericArrayField<ArrayType, D> &f, Vector<int, Dim> mesh, 
  std::function<datatype(datatype)> func)
      : varname(name), Field(f), Mesh(mesh), unitConvert(func) {}

  void writeBinary(const std::string &fName, int Overlap) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" " << "NumberOfComponents=\"" << D
      << "\">\n";
    f.close();

    // get submesh size
    std::size_t Size;
    Vector<int, Dim> SubMesh = Mesh - Vector<int, Dim>{2 * Overlap};
    if constexpr (Dim == 2) {
      Size = SubMesh[0] * SubMesh[1];
    } else if constexpr (Dim == 3) {
      Size = SubMesh[0] * SubMesh[1] * SubMesh[2];
    }
    std::size_t arrDSize = std::size_t(D * Size);

    // get data array for encoding
    datatype *data = new datatype[arrDSize];
    std::size_t id = 0;
    if constexpr (Dim == 2) {
      for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
        for (int i = Overlap; i < Mesh[0] - Overlap; ++i) {
          std::size_t arrId = i + j * Mesh[0];
          for (unsigned int vecid = 0; vecid < D; ++vecid) {
            data[id + vecid] = Field.get(arrId, vecid);
          }
          id += D;
        }
      }
    } else if constexpr (Dim == 3) {
      int XY = Mesh[0] * Mesh[1];
      for (int k = Overlap; k < Mesh[2] - Overlap; ++k) {
        for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
          for (int i = Overlap; i < Mesh[0] - Overlap; ++i) {
            std::size_t arrId = i + j * Mesh[0] + k * XY;
            for (unsigned int vecid = 0; vecid < D; ++vecid) {
              data[id + vecid] = Field.get(arrId, vecid);
            }
            id += D;
          }
        }
      }
    }

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(arrDSize * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);

    Base64Encoder<datatype> Encoder(fb, arrDSize);
    Encoder.encode(data, arrDSize);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

  void writeConvertedBinary(const std::string &fName, int Overlap) const override {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "format=\"binary\" encoding=\"base64\" " << "NumberOfComponents=\"" << D
      << "\">\n";
    f.close();

    // get submesh size
    std::size_t Size;
    Vector<int, Dim> SubMesh = Mesh - Vector<int, Dim>{2 * Overlap};
    if constexpr (Dim == 2) {
      Size = SubMesh[0] * SubMesh[1];
    } else if constexpr (Dim == 3) {
      Size = SubMesh[0] * SubMesh[1] * SubMesh[2];
    }
    std::size_t arrDSize = std::size_t(D * Size);

    // get data array for encoding
    datatype *data = new datatype[arrDSize];
    std::size_t id = 0;
    if constexpr (Dim == 2) {
      for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
        for (int i = Overlap; i < Mesh[0] - Overlap; ++i) {
          std::size_t arrId = i + j * Mesh[0];
          for (unsigned int vecid = 0; vecid < D; ++vecid) {
            data[id + vecid] = unitConvert(Field.get(arrId, vecid));
          }
          id += D;
        }
      }
    } else if constexpr (Dim == 3) {
      int XY = Mesh[0] * Mesh[1];
      for (int k = Overlap; k < Mesh[2] - Overlap; ++k) {
        for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
          for (int i = Overlap; i < Mesh[0] - Overlap; ++i) {
            std::size_t arrId = i + j * Mesh[0] + k * XY;
            for (unsigned int vecid = 0; vecid < D; ++vecid) {
              data[id + vecid] = unitConvert(Field.get(arrId, vecid));
            }
            id += D;
          }
        }
      }
    }

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    std::size_t binarySize = std::size_t(arrDSize * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);

    Base64Encoder<datatype> Encoder(fb, arrDSize);
    Encoder.encode(data, arrDSize);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

 private:
  std::string varname;
  // field data
  const GenericArrayField<ArrayType, D> &Field;
  // mesh info
  Vector<int, Dim> Mesh;
  // unit convert function pointer
  std::function<datatype(datatype)> unitConvert;
};

}  // namespace vtino