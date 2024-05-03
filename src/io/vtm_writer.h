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

// vtm_writer.h

#pragma once

#include "io/vti_writer.h"

// .pvd - .vtm - .vti

namespace vtmwriter {

// manager of generating a vtm file and series of vti files
template <typename T, unsigned int D>
class vtmWriter {
 private:
  std::string _dirname = "./vtkoutput/";
  std::string _vtidirname = "./vtkoutput/vtidata/";
  std::string _filename;
  std::vector<vtiwriter::vtiManager<T, D>> _vtiwriters;

  BlockGeometry<T, D>& BlockGeo;

 public:
  vtmWriter(std::string filename, BlockGeometry<T, D>& blockgeo)
      : _filename(filename), BlockGeo(blockgeo) {
    DirCreator::Create_Dir(_dirname);
    DirCreator::Create_Dir(_vtidirname);
    create_vtiwriters(blockgeo);
  }

  void Init() {
    _vtiwriters.clear();
    create_vtiwriters(BlockGeo);
  }

  void create_vtiwriters(BlockGeometry<T, D>& blockgeo) {
    for (const Block<T, D>& block : blockgeo.getBlocks()) {
      T voxsize = block.getVoxelSize();
      const Vector<int, D> ext = block.getMesh() - Vector<int, D>{1};
      const Vector<T, D> origin = block.getMinCenter();
      _vtiwriters.emplace_back(_filename, block.getBlockId(), voxsize, origin, ext);
    }
  }

  void addWriterSet(const vtiwriter::AbstWriterSet& writerset) {
    int i = 0;
    for (vtiwriter::vtiManager<T, D>& vti : _vtiwriters) {
      vti.addWriter(writerset.getWriter(i));
      ++i;
    }
  }
  template <typename... Args>
  void addWriterSet(const vtiwriter::AbstWriterSet& writerset, const Args&... args) {
    addWriterSet(writerset);
    addWriterSet(args...);
  }
  void Write() {
    // create vtm file
    std::string fullName = _vtidirname + _filename + ".vtm";
    vtmHeader(fullName);
    // create vti files
#pragma omp parallel for num_threads(Thread_Num)
    for (vtiwriter::vtiManager<T, D>& vti : _vtiwriters) {
      vti.Write();
      writevtm(fullName, vti.getFileName(), vti.getvtiBlockId());
    }
    // end vtm file
    vtmEnd(fullName);
  }
  void Write(int step) {
    std::string fullName = _vtidirname + _filename + std::to_string(step) + ".vtm";
    vtmHeader(fullName);
#pragma omp parallel for num_threads(Thread_Num)
    for (vtiwriter::vtiManager<T, D>& vti : _vtiwriters) {
      vti.Write(step);
      writevtm(fullName, vti.getFileName(step), vti.getvtiBlockId());
    }
    vtmEnd(fullName);
  }
  void WriteBinary() {
    mpi().barrier();
    // write vtis
#pragma omp parallel for num_threads(Thread_Num)
    for (vtiwriter::vtiManager<T, D>& vti : _vtiwriters) {
      vti.WriteBinary();
    }
    // write vtm
    std::string fullName = _vtidirname + _filename + ".vtm";
    IF_MPI_RANK(0) { vtmHeader(fullName); }
#ifdef MPI_ENABLED
    int vtisize = _vtiwriters.size();
    int totalvtisize;
    mpi().reduce(vtisize, totalvtisize, MPI_SUM);
    MPI_RANK(0)
    for (int i = 0; i < totalvtisize; ++i) {
      writevtm(fullName, getvtiFileName(i), i);
    }
#else
    for (int i = 0; i < _vtiwriters.size(); ++i) {
      writevtm(fullName, getvtiFileName(i), i);
    }
#endif
    vtmEnd(fullName);
  }

  void WriteBinary(int step) {
    mpi().barrier();
    // write vtis
#pragma omp parallel for num_threads(Thread_Num)
    for (vtiwriter::vtiManager<T, D>& vti : _vtiwriters) {
      vti.WriteBinary(step);
    }
    // write vtm
    std::string fullName = _vtidirname + _filename + std::to_string(step) + ".vtm";
    IF_MPI_RANK(0) { vtmHeader(fullName); }
    vtmHeader(fullName);
#ifdef MPI_ENABLED
    int vtisize = _vtiwriters.size();
    int totalvtisize;
    mpi().reduce(vtisize, totalvtisize, MPI_SUM);
    MPI_RANK(0)
    for (int i = 0; i < totalvtisize; ++i) {
      writevtm(fullName, getvtiFileName(i, step), i);
    }
#else
    for (int i = 0; i < _vtiwriters.size(); ++i) {
      writevtm(fullName, getvtiFileName(i, step), i);
    }
#endif
    vtmEnd(fullName);
  }

  std::string getvtiFileName(int id) {
    std::string fName = _filename + "_B" + std::to_string(id) + ".vti";
    return fName;
  }
  std::string getvtiFileName(int id, int step) {
    std::string fName =
      _filename + "_T" + std::to_string(step) + "_B" + std::to_string(id) + ".vti";
    return fName;
  }
  void vtmHeader(const std::string& fName) {
    std::ofstream f(fName, std::ios::trunc);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
      << "byte_order=\"LittleEndian\">\n"
      << "<vtkMultiBlockDataSet>\n";
    f.close();
  }
  void vtmEnd(const std::string& fName) {
    std::ofstream f(fName, std::ios::app);
    f << "</vtkMultiBlockDataSet>\n";
    f << "</VTKFile>\n";
    f.close();
  }
  void writevtm(const std::string& vtmName, const std::string& vtiName, int blockid) {
    std::ofstream f(vtmName, std::ios::app);
    f << "<Block index=\"" << blockid << "\" >\n";
    f << "<DataSet index= \"0\" "
      << "file=\"" << vtiName << "\">\n"
      << "</DataSet>\n";
    f << "</Block>\n";
    f.close();
  }
};

template <typename ArrayType>
class ScalerWriter : public vtiwriter::AbstWriterSet {
 private:
  std::vector<vtiwriter::ScalerWriter<ArrayType>> _scalerwriters;
  std::string VarName;

 public:
  template <typename FloatType, unsigned int Dim>
  ScalerWriter(
    std::string varname,
    const BlockFieldManager<GenericField<ArrayType, 1>, FloatType, Dim>& blockFM)
      : VarName(varname) {
    for (const BlockField<GenericField<ArrayType, 1>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _scalerwriters.emplace_back(varname, blockF.getField().getField(0));
    }
  }
  // create a writer for a single field
  ScalerWriter(std::string varname, GenericField<ArrayType, 1>& field)
      : VarName(varname) {
    _scalerwriters.emplace_back(varname, field.getField(0));
  }

  template <typename FloatType, unsigned int Dim>
  void Init(
    const BlockFieldManager<GenericField<ArrayType, 1>, FloatType, Dim>& blockFM) {
    _scalerwriters.clear();
    for (const BlockField<GenericField<ArrayType, 1>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _scalerwriters.emplace_back(VarName, blockF.getField().getField(0));
    }
  }

  const vtiwriter::AbstractWriter& getWriter(int i) const override {
    return _scalerwriters[i];
  }
};

// GenericArray<Vector<datatype, D>> &Array
template <typename datatype, unsigned int D>
class VectorWriter : public vtiwriter::AbstWriterSet {
 private:
  std::vector<vtiwriter::VectorWriter<datatype, D>> _vectorwriters;
  std::string VarName;

 public:
  template <typename FloatType, unsigned int Dim>
  VectorWriter(
    std::string varname,
    const BlockFieldManager<VectorFieldAOS<datatype, D>, FloatType, Dim>& blockFM)
      : VarName(varname) {
    for (const BlockField<VectorFieldAOS<datatype, D>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _vectorwriters.emplace_back(varname, blockF.getField().getField(0));
    }
  }
  VectorWriter(std::string varname, std::vector<VectorFieldAOS<datatype, D>*> field)
      : VarName(varname) {
    for (int i = 0; i < field.size(); ++i) {
      _vectorwriters.emplace_back(varname, field[i]->getField(0));
    }
  }
  // create a writer for a single field
  VectorWriter(std::string varname, VectorFieldAOS<datatype, D>& field)
      : VarName(varname) {
    _vectorwriters.emplace_back(varname, field.getField(0));
  }

  template <typename FloatType, unsigned int Dim>
  void Init(
    const BlockFieldManager<VectorFieldAOS<datatype, D>, FloatType, Dim>& blockFM) {
    _vectorwriters.clear();
    for (const BlockField<VectorFieldAOS<datatype, D>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _vectorwriters.emplace_back(VarName, blockF.getField().getField(0));
    }
  }

  const vtiwriter::AbstractWriter& getWriter(int i) const override {
    return _vectorwriters[i];
  }
};

// GenericField<ArrayType, D>
template <typename ArrayType, unsigned int D>
class VectorSOAWriter : public vtiwriter::AbstWriterSet {
 private:
  std::vector<vtiwriter::VectorSOAWriter<ArrayType, D>> _vectorsoawriters;
  std::string VarName;

 public:
  template <typename FloatType, unsigned int Dim>
  VectorSOAWriter(
    std::string varname,
    const BlockFieldManager<GenericField<ArrayType, D>, FloatType, Dim>& blockFM)
      : VarName(varname) {
    for (const BlockField<GenericField<ArrayType, D>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _vectorsoawriters.emplace_back(varname, blockF.getField());
    }
  }
  // create a writer for a single field
  VectorSOAWriter(std::string varname, GenericField<ArrayType, D>& field)
      : VarName(varname) {
    _vectorsoawriters.emplace_back(varname, &field);
  }

  template <typename FloatType, unsigned int Dim>
  void Init(
    const BlockFieldManager<GenericField<ArrayType, D>, FloatType, Dim>& blockFM) {
    _vectorsoawriters.clear();
    for (const BlockField<GenericField<ArrayType, D>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _vectorsoawriters.emplace_back(VarName, blockF.getField());
    }
  }

  const vtiwriter::AbstractWriter& getWriter(int i) const override {
    return _vectorsoawriters[i];
  }
};

}  // namespace vtmwriter

// vtm writer with handling overlap
namespace vtmo {

// manager of generating a vtm file and series of vti files
template <typename T, unsigned int D>
class vtmWriter {
 private:
  std::string _dirname = "./vtkoutput/";
  std::string _vtidirname = "./vtkoutput/vtidata/";
  std::string _filename;
  std::vector<vtino::vtiManager<T, D>> _vtiwriters;
  // if block overlap is larger than this value, use overlap - threshold
  // if block overlap is smaller than or equal to this value, use 0
  // you could set it as:
  // default: -1, use default overlap, vtmo::vtmWriter becomes vtmwriter::vtmWriter
  // 0, no overlap, not recommended for paraview-surface rendering
  // 1, each block has 1 overlap
  int _Overlap_Threshold;
  BlockGeometry<T, D>& BlockGeo;

 public:
  vtmWriter(std::string filename, BlockGeometry<T, D>& blockgeo, int overlapth = -1)
      : _filename(filename), _Overlap_Threshold(overlapth), BlockGeo(blockgeo) {
    DirCreator::Create_Dir(_dirname);
    DirCreator::Create_Dir(_vtidirname);
    create_vtiwriters(blockgeo);
  }

  void Init() {
    _vtiwriters.clear();
    create_vtiwriters(BlockGeo);
  }

  void create_vtiwriters(BlockGeometry<T, D>& blockgeo) {
    for (const Block<T, D>& block : blockgeo.getBlocks()) {
      T voxsize = block.getVoxelSize();
      const Vector<int, D> ext = block.getMesh() - Vector<int, D>{1};
      const Vector<T, D> origin = block.getMinCenter();
      int overlap;
      if (_Overlap_Threshold == -1) {
        overlap = 0;
      } else {
        if (_Overlap_Threshold < block.getOverlap()) {
          overlap = block.getOverlap() - _Overlap_Threshold;
        } else {
          overlap = 0;
        }
      }
      _vtiwriters.emplace_back(_filename, block.getBlockId(), voxsize, origin, ext, overlap);
    }
  }

  void addWriterSet(const vtino::AbstWriterSet& writerset) {
    int i = 0;
    for (vtino::vtiManager<T, D>& vti : _vtiwriters) {
      vti.addWriter(writerset.getWriter(i));
      ++i;
    }
  }
  template <typename... Args>
  void addWriterSet(const vtino::AbstWriterSet& writerset, const Args&... args) {
    addWriterSet(writerset);
    addWriterSet(args...);
  }

  void WriteBinary() {
    mpi().barrier();
    // write vtis
#pragma omp parallel for num_threads(Thread_Num)
    for (vtino::vtiManager<T, D>& vti : _vtiwriters) {
      vti.WriteBinary();
    }
    // write vtm
    std::string fullName = _vtidirname + _filename + ".vtm";
    IF_MPI_RANK(0) { vtmHeader(fullName); }
#ifdef MPI_ENABLED
    int vtisize = _vtiwriters.size();
    int totalvtisize;
    mpi().reduce(vtisize, totalvtisize, MPI_SUM);
    MPI_RANK(0)
    for (int i = 0; i < totalvtisize; ++i) {
      writevtm(fullName, getvtiFileName(i), i);
    }
#else
    for (int i = 0; i < _vtiwriters.size(); ++i) {
      writevtm(fullName, getvtiFileName(i), i);
    }
#endif
    vtmEnd(fullName);
  }

  void WriteBinary(int step) {
    mpi().barrier();
    // write vtis
#pragma omp parallel for num_threads(Thread_Num)
    for (vtino::vtiManager<T, D>& vti : _vtiwriters) {
      vti.WriteBinary(step);
    }
    // write vtm
    std::string fullName = _vtidirname + _filename + std::to_string(step) + ".vtm";
    IF_MPI_RANK(0) { vtmHeader(fullName); }
    vtmHeader(fullName);
#ifdef MPI_ENABLED
    int vtisize = _vtiwriters.size();
    int totalvtisize;
    mpi().reduce(vtisize, totalvtisize, MPI_SUM);
    MPI_RANK(0)
    for (int i = 0; i < totalvtisize; ++i) {
      writevtm(fullName, getvtiFileName(i, step), i);
    }
#else
    for (int i = 0; i < _vtiwriters.size(); ++i) {
      writevtm(fullName, getvtiFileName(i, step), i);
    }
#endif
    vtmEnd(fullName);
  }

  std::string getvtiFileName(int id) {
    std::string fName = _filename + "_B" + std::to_string(id) + ".vti";
    return fName;
  }
  std::string getvtiFileName(int id, int step) {
    std::string fName =
      _filename + "_T" + std::to_string(step) + "_B" + std::to_string(id) + ".vti";
    return fName;
  }
  void vtmHeader(const std::string& fName) {
    std::ofstream f(fName, std::ios::trunc);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
      << "byte_order=\"LittleEndian\">\n"
      << "<vtkMultiBlockDataSet>\n";
    f.close();
  }
  void vtmEnd(const std::string& fName) {
    std::ofstream f(fName, std::ios::app);
    f << "</vtkMultiBlockDataSet>\n";
    f << "</VTKFile>\n";
    f.close();
  }
  void writevtm(const std::string& vtmName, const std::string& vtiName, int blockid) {
    std::ofstream f(vtmName, std::ios::app);
    f << "<Block index=\"" << blockid << "\" >\n";
    f << "<DataSet index= \"0\" "
      << "file=\"" << vtiName << "\">\n"
      << "</DataSet>\n";
    f << "</Block>\n";
    f.close();
  }
};

template <typename ArrayType, unsigned int Dim>
class ScalerWriter : public vtino::AbstWriterSet {
 private:
  std::vector<vtino::ScalerWriter<ArrayType, Dim>> _scalerwriters;
  std::string VarName;

 public:
  ScalerWriter(std::string varname, std::vector<GenericField<ArrayType, 1>*> field,
               std::vector<Vector<int, Dim>*> meshes)
      : VarName(varname) {
    for (int i = 0; i < field.size(); ++i) {
      _scalerwriters.emplace_back(varname, field[i]->getField(0), *meshes[i]);
    }
  }
  template <typename FloatType>
  ScalerWriter(
    std::string varname,
    const BlockFieldManager<GenericField<ArrayType, 1>, FloatType, Dim>& blockFM)
      : VarName(varname) {
    for (const BlockField<GenericField<ArrayType, 1>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _scalerwriters.emplace_back(varname, blockF.getField().getField(0),
                                  blockF.getBlock().getMesh());
    }
  }

  template <typename FloatType>
  void Init(
    const BlockFieldManager<GenericField<ArrayType, 1>, FloatType, Dim>& blockFM) {
    _scalerwriters.clear();
    for (const BlockField<GenericField<ArrayType, 1>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _scalerwriters.emplace_back(VarName, blockF.getField().getField(0),
                                  blockF.getBlock().getMesh());
    }
  }

  const vtino::AbstractWriter& getWriter(int i) const override {
    return _scalerwriters[i];
  }
};

template <typename datatype, unsigned int Dim, unsigned int D>
class VectorWriter : public vtino::AbstWriterSet {
 private:
  std::vector<vtino::VectorWriter<datatype, Dim, D>> _vectorwriters;
  std::string VarName;

 public:
  VectorWriter(std::string varname, std::vector<VectorFieldAOS<datatype, D>*> field,
               std::vector<Vector<int, Dim>*> meshes)
      : VarName(varname) {
    for (int i = 0; i < field.size(); ++i) {
      _vectorwriters.emplace_back(varname, field[i]->getField(0), *meshes[i]);
    }
  }
  template <typename FloatType>
  VectorWriter(
    std::string varname,
    const BlockFieldManager<VectorFieldAOS<datatype, D>, FloatType, Dim>& blockFM)
      : VarName(varname) {
    for (const BlockField<VectorFieldAOS<datatype, D>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _vectorwriters.emplace_back(varname, blockF.getField().getField(0),
                                  blockF.getBlock().getMesh());
    }
  }

  template <typename FloatType>
  void Init(
    const BlockFieldManager<VectorFieldAOS<datatype, D>, FloatType, Dim>& blockFM) {
    _vectorwriters.clear();
    for (const BlockField<VectorFieldAOS<datatype, D>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _vectorwriters.emplace_back(VarName, blockF.getField().getField(0),
                                  blockF.getBlock().getMesh());
    }
  }

  const vtino::AbstractWriter& getWriter(int i) const override {
    return _vectorwriters[i];
  }
};

template <typename ArrayType, unsigned int Dim, unsigned int D>
class VectorSOAWriter : public vtino::AbstWriterSet {
 private:
  std::vector<vtino::VectorSOAWriter<ArrayType, Dim, D>> _vectorsoawriters;
  std::string VarName;

 public:
  VectorSOAWriter(std::string varname, std::vector<GenericField<ArrayType, D>*> field,
                  std::vector<Vector<int, Dim>*> meshes)
      : VarName(varname) {
    for (int i = 0; i < field.size(); ++i) {
      _vectorsoawriters.emplace_back(varname, &field[i], *meshes[i]);
    }
  }
  template <typename FloatType>
  VectorSOAWriter(
    std::string varname,
    const BlockFieldManager<GenericField<ArrayType, D>, FloatType, Dim>& blockFM)
      : VarName(varname) {
    for (const BlockField<GenericField<ArrayType, D>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _vectorsoawriters.emplace_back(varname, blockF.getField(),
                                     blockF.getBlock().getMesh());
    }
  }

  template <typename FloatType>
  void Init(
    const BlockFieldManager<GenericField<ArrayType, D>, FloatType, Dim>& blockFM) {
    _vectorsoawriters.clear();
    for (const BlockField<GenericField<ArrayType, D>, FloatType, Dim>& blockF :
         blockFM.getBlockFields()) {
      _vectorsoawriters.emplace_back(VarName, blockF.getField(),
                                     blockF.getBlock().getMesh());
    }
  }

  const vtino::AbstractWriter& getWriter(int i) const override {
    return _vectorsoawriters[i];
  }
};

}  // namespace vtmo