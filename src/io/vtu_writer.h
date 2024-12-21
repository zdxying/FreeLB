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

// vtu_writer.h

#pragma once

// store any callable object and used in unit converter when writing vti files
#include <functional>

#include "data_struct/Vector.h"
#include "data_struct/block_lattice.h"
#include "io/base64.h"
#include "io/basic_writer.h"

// offlattice
#include "offLattice/triangleSet.h"

// https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#unstructuredgrid
/*
<VTKFile type="UnstructuredGrid" ...>
  <UnstructuredGrid>
    <Piece NumberOfPoints="#" NumberOfCells="#">
    <PointData>...</PointData>
    <CellData>...</CellData>
    <Points>...</Points>
    <Cells>...</Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
*/

namespace vtuSurface {

class AbstractWriter {
  public:
    virtual void writeBinary(const std::string &fName) const = 0;
    virtual void write(const std::string &fName) const = 0;
};


// manager of generating a single vtu file
template <typename T>
class vtuManager {
 private:
  std::string _dirname = "./vtkoutput/vtudata/";
  std::string _filename;

  const offlat::TriangleSet<T>& _triangleSet;

  std::vector<const AbstractWriter *> _writers;

  void WriteBinaryImpl(const std::string &fName) {
    mpi().barrier();
    vtuHeader(fName);
    vtuPointsBinary(fName);
    vtuCellsBinary(fName);
    for (const AbstractWriter *writer : _writers) {
      writer->writeBinary(fName);
    }
    vtuEnd(fName);
  }

  void WriteImpl(const std::string &fName) {
    mpi().barrier();
    vtuHeader(fName);
    vtuPoints(fName);
    vtuCells(fName);
    for (const AbstractWriter *writer : _writers) {
      writer->write(fName);
    }
    vtuEnd(fName);
  }

  void WritevtmImpl() {
#ifdef MPI_ENABLED
    IF_MPI_RANK(0) { 
      std::string fullName = _dirname + _filename + ".vtm";
      vtmHeader(fullName);     
      for (int iRank = 0; iRank < mpi().getSize(); ++iRank) {
        vtmContent(fullName, vtuFileName(iRank), iRank);
      }
      vtmEnd(fullName); 
    }
#endif
  }
  void WritevtmImpl(int step) {
#ifdef MPI_ENABLED
    IF_MPI_RANK(0) { 
      std::string fullName = _dirname + _filename + "_T" + std::to_string(step) + ".vtm";
      vtmHeader(fullName);     
      for (int iRank = 0; iRank < mpi().getSize(); ++iRank) {
        vtmContent(fullName, vtuFileName(iRank, step), iRank);
      }
      vtmEnd(fullName); 
    }
#endif
  }

  void vtuPointsBinary(const std::string &fName) {
    const std::size_t _numPoints = _triangleSet.getTriangles().size() * 3;
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<T>(type);
    
    f << "<Points>\n";
    f << "<DataArray type=\"" << type <<"\" NumberOfComponents=\"" << 3 
      << "\" format=\"binary\" " << "encoding=\"base64\">\n";
    f.close();

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    unsigned int uintBinarySize = static_cast<unsigned int>(_numPoints * 3 * sizeof(T));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<T> Encoder(fb, _numPoints * 3);    
    for (const auto &triangle : _triangleSet.getTriangles()) {
      for (unsigned int i = 0; i < 3; ++i) {
        Encoder.encode(triangle[i].data(), 3);
      }
    }
    fb.close();

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff << "</Points>\n";
    ff.close();
  }

  void vtuCellsBinary(const std::string &fName) {
    const std::size_t _numCells = _triangleSet.getTriangles().size();
    std::ofstream f(fName, std::ios::out | std::ios::app);
    f << "<Cells>\n";

    // ---------connectivity

    f << "<DataArray type=\"UInt32\" Name=\"connectivity\" "
      << "format=\"binary\" encoding=\"base64\">\n";
    f.close();

    // generate an array of connectivity: 0, 1, 2, 3, 4, ...
    // size = _numCells * 3
    unsigned int* arr = new unsigned int[_numCells * 3];
    for (std::size_t i = 0; i < _numCells * 3; ++i) {
      arr[i] = static_cast<unsigned int>(i);
    }

    std::ofstream fcb(fName, std::ios::out | std::ios::app | std::ios::binary);
    unsigned int uintBinarySize = static_cast<unsigned int>(_numCells * 3 * sizeof(unsigned int));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoderc(fcb, 1);
    sizeEncoderc.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<unsigned int> Encoderc(fcb, _numCells * 3);
    Encoderc.encode(arr, _numCells * 3);
    fcb.close();
    
    
    std::ofstream fo(fName, std::ios::out | std::ios::app);
    fo << "\n</DataArray>\n";

    // ---------offsets

    fo << "<DataArray type=\"UInt32\" Name=\"offsets\" "
       << "format=\"binary\" encoding=\"base64\">\n";
    fo.close();

    // generate an array of offsets: 3, 6, 9, 12, 15, ...
    // size = _numCells, use existing array arr
    for (std::size_t i = 0; i < _numCells; ++i) {
      arr[i] = static_cast<unsigned int>((i + 1) * 3);
    }

    std::ofstream fob(fName, std::ios::out | std::ios::app | std::ios::binary);
    uintBinarySize = static_cast<unsigned int>(_numCells * sizeof(unsigned int));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncodero(fob, 1);
    sizeEncodero.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<unsigned int> Encoderro(fob, _numCells);
    Encoderro.encode(arr, _numCells);
    fob.close();
    delete[] arr;

    std::ofstream ft(fName, std::ios::out | std::ios::app);
    ft << "\n</DataArray>\n";

    // ---------types
    
    ft << "<DataArray type=\"UInt8\" Name=\"types\" "
       << "format=\"binary\" encoding=\"base64\">\n";
    ft.close();

    // generate an array of type VTK_TRIANGLE: 5, 5, 5, 5, 5, ...
    // size = _numCells
    std::uint8_t* arrt = new std::uint8_t[_numCells];
    std::fill(arrt, arrt + _numCells, std::uint8_t{5});

    std::ofstream ftb(fName, std::ios::out | std::ios::app | std::ios::binary);
    uintBinarySize = static_cast<unsigned int>(_numCells * sizeof(std::uint8_t));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncodert(ftb, 1);
    sizeEncodert.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<std::uint8_t> Encodert(ftb, _numCells);
    Encodert.encode(arrt, _numCells);
    ftb.close();
    delete[] arrt;

    std::ofstream fft(fName, std::ios::out | std::ios::app);
    fft << "\n</DataArray>\n";
    fft << "</Cells>\n";

    // ---------point data

    fft << "<PointData Scalars=\"scalars\" Vectors=\"vectors\">\n";
    fft.close();
  }

  void vtuPoints(const std::string &fName) {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<T>(type);
    
    f << "<Points>\n";
    f << "<DataArray type=\"" << type <<"\" NumberOfComponents=\"" << 3 
      << "\" format=\"ascii\">\n";
    for (const auto &triangle : _triangleSet.getTriangles()) {
      for (unsigned int i = 0; i < 3; ++i) {
        f << triangle[i][0] << " " << triangle[i][1] << " " << triangle[i][2] << " ";
      }
    }
    f << "\n</DataArray>\n";
    f << "</Points>\n";
    f.close();
  }

  void vtuCells(const std::string &fName) {
    const std::size_t _numCells = _triangleSet.getTriangles().size();
    std::ofstream f(fName, std::ios::out | std::ios::app);
    f << "<Cells>\n";

    // ---------connectivity
    f << "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";
    // generate an array of connectivity: 0, 1, 2, 3, 4, ...
    for (unsigned int i = 0; i < _numCells; ++i) {
      f << i * 3 << " " << i * 3 + 1 << " " << i * 3 + 2 << " ";
    }
    f << "\n</DataArray>\n";

    // ---------offsets
    f << "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n";
    // generate an array of offsets: 3, 6, 9, 12, 15, ...
    for (unsigned int i = 0; i < _numCells; ++i) {
      f << (i + 1) * 3 << " ";
    }
    f << "\n</DataArray>\n";

    // ---------types
    f << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    // generate an array of type VTK_TRIANGLE: 5, 5, 5, 5, 5, ...
    for (std::size_t i = 0; i < _numCells; ++i) {
      f << 5 << " ";
    }
    f << "\n</DataArray>\n";
    f << "</Cells>\n";

    // ---------point data
    f << "<PointData Scalars=\"scalars\" Vectors=\"vectors\">\n";
    f.close();
  }

  void vtuHeader(const std::string &fName) {
    const std::size_t _numCells = _triangleSet.getTriangles().size();
    const std::size_t _numPoints = _triangleSet.getTriangles().size() * 3;
    std::ofstream f(fName);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
      << "byte_order=\"LittleEndian\">\n";
    f << "<UnstructuredGrid>\n";
    f << "<Piece NumberOfPoints=\"" << _numPoints << "\" NumberOfCells=\"" << _numCells << "\">\n";
    f.close();
  }
  void vtuEnd(const std::string &fName) {
    std::ofstream f(fName, std::ios::app);
    f << "</PointData>\n";
    f << "</Piece>\n";
    f << "</UnstructuredGrid>\n";
    f << "</VTKFile>\n";
    f.close();
  }

  void vtmHeader(const std::string& fName) {
    std::ofstream f(fName, std::ios::out | std::ios::trunc);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
      << "byte_order=\"LittleEndian\">\n"
      << "<vtkMultiBlockDataSet>\n";
    f.close();
  }
  void vtmEnd(const std::string& fName) {
    std::ofstream f(fName, std::ios::out | std::ios::app);
    f << "</vtkMultiBlockDataSet>\n";
    f << "</VTKFile>\n";
    f.close();
  }
  void vtmContent(const std::string& vtmName, const std::string& vtuName, int rank) {
    std::ofstream f(vtmName, std::ios::out | std::ios::app);
    f << "<Block index=\"" << rank << "\" >\n";
    f << "<DataSet index= \"0\" " << "file=\"" << vtuName << "\">\n"
      << "</DataSet>\n";
    f << "</Block>\n";
    f.close();
  }

  std::string vtuFileName(int rank) {
    std::string fName{};
#ifdef MPI_ENABLED
    fName = _filename + "_B" + std::to_string(rank) + ".vtu";
#else
    fName = _filename + ".vtu";
#endif
    return fName;
  }
  std::string vtuFileName(int rank, int step) {
    std::string fName{};
#ifdef MPI_ENABLED
    fName = _filename + "_T" + std::to_string(step) 
    + "_B" + std::to_string(rank)+ ".vtu";
#else
    fName = _filename + "_T" + std::to_string(step) + ".vtu";
#endif
    return fName;
  }
  std::string vtuFullName() {
    std::string fName{};
#ifdef MPI_ENABLED
    fName = _dirname + _filename + "_B" + std::to_string(mpi().getRank()) + ".vtu";
#else
    fName = _dirname + _filename + ".vtu";
#endif
    return fName;
  }
  std::string vtuFullName(int step) {
    std::string fName{};
#ifdef MPI_ENABLED
    fName = _dirname + _filename + "_T" + std::to_string(step) 
    + "_B" + std::to_string(mpi().getRank())+ ".vtu";
#else
    fName = _dirname + _filename + "_T" + std::to_string(step) + ".vtu";
#endif
    return fName;
  }

 public:
  vtuManager(const std::string& filename, const offlat::TriangleSet<T>& triangleSet)
      : _filename(filename), _triangleSet(triangleSet) {
    DirCreator::Create_Dir("./vtkoutput/");
    DirCreator::Create_Dir(_dirname);
  }

  void addWriter(const AbstractWriter &writer) { _writers.push_back(&writer); }
  template <typename... Args>
  void addWriter(const AbstractWriter &writer, const Args &...args) {
    addWriter(writer);
    addWriter(args...);
  }

  void WriteBinary() {
    WriteBinaryImpl(vtuFullName());
    WritevtmImpl();
  }
  void WriteBinary(int step) {
    WriteBinaryImpl(vtuFullName(step));
    WritevtmImpl(step);
  }

  void Write() {
    WriteImpl(vtuFullName());
    WritevtmImpl();
  }
  void Write(int step) {
    WriteImpl(vtuFullName(step));
    WritevtmImpl(step);
  }
};


// get filed data from field manager, filterd by flag manager
template <typename FieldType, typename FloatType>
void getFieldData(const BlockFieldManager<FieldType, FloatType, 3>& fieldManager, 
  const std::vector<std::vector<offlat::TriangleIdx<FloatType>>>& triIdxvecs, 
  typename FieldType::value_type* data) {
  
  std::size_t id{};
  for (std::size_t iblock = 0; iblock < fieldManager.size(); ++iblock) {
    const auto& field = fieldManager.getBlockField(iblock).getField(0);
    const std::vector<offlat::TriangleIdx<FloatType>>& triIdxs = triIdxvecs[iblock];

    for (const offlat::TriangleIdx<FloatType>& triIdx : triIdxs) {
      data[id] = triIdx[0].getIntp(field[triIdx[0].id1], field[triIdx[0].id2]);
      ++id;
      data[id] = triIdx[1].getIntp(field[triIdx[1].id1], field[triIdx[1].id2]);
      ++id;
      data[id] = triIdx[2].getIntp(field[triIdx[2].id1], field[triIdx[2].id2]);
      ++id;
    }

  }
}

// get filed data from field manager, filterd by flag manager, convert by unit converter
template <typename FieldType, typename FloatType>
void getFieldData(const BlockFieldManager<FieldType, FloatType, 3>& fieldManager, 
  const std::vector<std::vector<offlat::TriangleIdx<FloatType>>>& triIdxvecs, 
  typename FieldType::value_type* data,
  std::function<typename FieldType::value_type(typename FieldType::value_type)> f) {
  
  std::size_t id{};
  for (std::size_t iblock = 0; iblock < fieldManager.size(); ++iblock) {
    const auto& field = fieldManager.getBlockField(iblock).getField(0);
    const std::vector<offlat::TriangleIdx<FloatType>>& triIdxs = triIdxvecs[iblock];

    for (const offlat::TriangleIdx<FloatType>& triIdx : triIdxs) {
      data[id] = f(triIdx[0].getIntp(field[triIdx[0].id1], field[triIdx[0].id2]));
      ++id;
      data[id] = f(triIdx[1].getIntp(field[triIdx[1].id1], field[triIdx[1].id2]));
      ++id;
      data[id] = f(triIdx[2].getIntp(field[triIdx[2].id1], field[triIdx[2].id2]));
      ++id;
    }

  }
}

template <typename FieldType, typename FloatType>
class ScalarWriter : public AbstractWriter {
 public:
  using datatype = typename FieldType::value_type;

  ScalarWriter(const std::string& name, const BlockFieldManager<FieldType, FloatType, 3> &fieldManager, 
    const offlat::TriangleSet<FloatType>& triSet)
      : varname(name), _fieldManager(fieldManager), _triSet(triSet) {}

  void writeBinary(const std::string &fName) const override {
    // number of point data
    std::size_t _pointNum = _triSet.getTriangles().size()*3;
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"1\" "
      << "format=\"binary\" encoding=\"base64\">\n";
    f.close();

    datatype *data = new datatype[_pointNum];
    getFieldData<FieldType, FloatType>(_fieldManager, _triSet.getTriangleIdxs(), data);

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    unsigned int uintBinarySize = static_cast<unsigned int>(_pointNum * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<datatype> Encoder(fb, _pointNum);
    Encoder.encode(data, _pointNum);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

  void write(const std::string &fName) const override {
    // number of point data
    std::size_t _pointNum = _triSet.getTriangles().size()*3;
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"1\" format=\"ascii\">\n";

    datatype *data = new datatype[_pointNum];
    getFieldData<FieldType, FloatType>(_fieldManager, _triSet.getTriangleIdxs(), data);

    for (std::size_t i = 0; i < _pointNum; ++i) {
      f << data[i] << " ";
    }

    f << "\n</DataArray>\n";
    f.close();
    delete[] data;
  }

 private:
  std::string varname;
  // field data
  const BlockFieldManager<FieldType, FloatType, 3>& _fieldManager;
  // triangle info
  const offlat::TriangleSet<FloatType>& _triSet;
};

template <typename FieldType, typename FloatType>
class physScalarWriter : public AbstractWriter {
 public:
  using datatype = typename FieldType::value_type;

  // std::bind(&uintConvclass::func, &unitConv, std::placeholders::_1); or 
  // [&unitConv](T x) { return unitConv.func(x); };
  physScalarWriter(const std::string& name, const BlockFieldManager<FieldType, FloatType, 3> &fieldManager, 
    const offlat::TriangleSet<FloatType>& triSet, std::function<datatype(datatype)> func)
      : varname(name), _fieldManager(fieldManager), _triSet(triSet), unitConvert(func) {}

  void writeBinary(const std::string &fName) const override {
    // number of point data
    std::size_t _pointNum = _triSet.getTriangles().size()*3;    
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"1\" "
      << "format=\"binary\" encoding=\"base64\">\n";
    f.close();

    datatype *data = new datatype[_pointNum];
    getFieldData<FieldType, FloatType>(_fieldManager, _triSet.getTriangleIdxs(), data, unitConvert);

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    unsigned int uintBinarySize = static_cast<unsigned int>(_pointNum * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<datatype> Encoder(fb, _pointNum);
    Encoder.encode(data, _pointNum);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

  void write(const std::string &fName) const override {
    // number of point data
    std::size_t _pointNum = _triSet.getTriangles().size()*3;    
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"1\" format=\"ascii\">\n";

    datatype *data = new datatype[_pointNum];
    getFieldData<FieldType, FloatType>(_fieldManager, _triSet.getTriangleIdxs(), data, unitConvert);

    for (std::size_t i = 0; i < _pointNum; ++i) {
      f << data[i] << " ";
    }

    f << "\n</DataArray>\n";
    f.close();
    delete[] data;
  }

 private:
  std::string varname;
  // field data
  const BlockFieldManager<FieldType, FloatType, 3>& _fieldManager;
  // triangle info
  const offlat::TriangleSet<FloatType>& _triSet;
  // unit convert function pointer
  std::function<datatype(datatype)> unitConvert;
};


template <typename FieldType, typename FloatType>
class VectorWriter : public AbstractWriter {
 public:
  using vectortype = typename FieldType::value_type;
  using datatype = typename vectortype::value_type;
  static constexpr unsigned int D = vectortype::vector_dim;

  VectorWriter(const std::string& name, const BlockFieldManager<FieldType, FloatType, 3> &fieldManager, 
    const offlat::TriangleSet<FloatType>& triSet)
      : varname(name), _fieldManager(fieldManager), _triSet(triSet) {}

  void writeBinary(const std::string &fName) const override {
    // number of point data
    std::size_t _pointNum = _triSet.getTriangles().size()*3;
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << D << "\" "
      << "format=\"binary\" encoding=\"base64\">\n";
    f.close();

    vectortype *data = new vectortype[_pointNum];
    getFieldData<FieldType, FloatType>(_fieldManager, _triSet.getTriangleIdxs(), data);

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    unsigned int uintBinarySize = static_cast<unsigned int>(_pointNum * D * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<datatype> Encoder(fb, _pointNum * D);
    Encoder.encode(data->data(), _pointNum * D);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

  void write(const std::string &fName) const override {
    // number of point data
    std::size_t _pointNum = _triSet.getTriangles().size()*3;
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << D <<"\" format=\"ascii\">\n";

    vectortype *data = new vectortype[_pointNum];
    getFieldData<FieldType, FloatType>(_fieldManager, _triSet.getTriangleIdxs(), data);

    for (std::size_t i = 0; i < _pointNum; ++i) {
      for (unsigned int j = 0; j < D; ++j) {
        f << data[i][j] << " ";
      }
    }

    f << "\n</DataArray>\n";
    f.close();
    delete[] data;
  }

 private:
  std::string varname;
  // field data
  const BlockFieldManager<FieldType, FloatType, 3>& _fieldManager;
  // triangle info
  const offlat::TriangleSet<FloatType>& _triSet;
};

template <typename FieldType, typename FloatType>
class physVectorWriter : public AbstractWriter {
 public:
  using vectortype = typename FieldType::value_type;
  using datatype = typename vectortype::value_type;
  static constexpr unsigned int D = vectortype::vector_dim;

  // std::bind(&uintConvclass::func, &unitConv, std::placeholders::_1); or 
  // [&unitConv](T x) { return unitConv.func(x); };
  physVectorWriter(const std::string& name, const BlockFieldManager<FieldType, FloatType, 3> &fieldManager, 
    const offlat::TriangleSet<FloatType>& triSet, std::function<vectortype(vectortype)> func)
      : varname(name), _fieldManager(fieldManager), _triSet(triSet), unitConvert(func) {}

  void writeBinary(const std::string &fName) const override {
    // number of point data
    std::size_t _pointNum = _triSet.getTriangles().size()*3;    
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << D << "\" "
      << "format=\"binary\" encoding=\"base64\">\n";
    f.close();

    vectortype *data = new vectortype[_pointNum];
    getFieldData<FieldType, FloatType>(_fieldManager, _triSet.getTriangleIdxs(), data, unitConvert);

    std::ofstream fb(fName, std::ios::out | std::ios::app | std::ios::binary);
    unsigned int uintBinarySize = static_cast<unsigned int>(_pointNum * D * sizeof(datatype));
    // writes first number, the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(fb, 1);
    sizeEncoder.encode(&uintBinarySize, 1);
    // writes the data
    Base64Encoder<datatype> Encoder(fb, _pointNum * D);
    Encoder.encode(data->data(), _pointNum * D);
    fb.close();
    delete[] data;

    std::ofstream ff(fName, std::ios::out | std::ios::app);
    ff << "\n</DataArray>\n";
    ff.close();
  }

  void write(const std::string &fName) const override {
    // number of point data
    std::size_t _pointNum = _triSet.getTriangles().size()*3;    
    std::ofstream f(fName, std::ios::out | std::ios::app);
    std::string type;
    getVTKTypeString<datatype>(type);
    f << "<DataArray type=\"" << type << "\" Name=\"" << varname << "\" "
      << "NumberOfComponents=\"" << D <<"\" format=\"ascii\">\n";

    vectortype *data = new vectortype[_pointNum];
    getFieldData<FieldType, FloatType>(_fieldManager, _triSet.getTriangleIdxs(), data, unitConvert);

    for (std::size_t i = 0; i < _pointNum; ++i) {
      for (unsigned int j = 0; j < D; ++j) {
        f << data[i][j] << " ";
      }
    }

    f << "\n</DataArray>\n";
    f.close();
    delete[] data;
  }

 private:
  std::string varname;
  // field data
  const BlockFieldManager<FieldType, FloatType, 3>& _fieldManager;
  // triangle info
  const offlat::TriangleSet<FloatType>& _triSet;
  // unit convert function pointer
  std::function<datatype(datatype)> unitConvert;
};

}  // namespace vtuSurface