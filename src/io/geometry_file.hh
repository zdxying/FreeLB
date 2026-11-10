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

// geometry_file.hh

#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// AABB::check() uses MPI_RANK; include the macro header before the geometry
// implementation, as io/voxel_geometry_builder.hh does for StlReader.
#include "utils/mpi_macros.h"
#include "geometry/voxel_geometry.hh"
#include "io/geometry_file.h"

namespace flbgeo {
namespace detail {

constexpr std::uint32_t kGeometryFileVersion = 1;
constexpr std::size_t kMaxPathLength = 4096;
constexpr char kGeometryFileMagic[8] = {'F', 'L', 'B', 'G', 'E', 'O', '\0', '\0'};
constexpr std::uint64_t kFnvOffset = 14695981039346656037ull;
constexpr std::uint64_t kFnvPrime = 1099511628211ull;

[[noreturn]] inline void fail(const std::string& message) {
  std::cerr << "[flbgeo] geometry file: " << message << std::endl;
  std::exit(1);
}

inline std::uint64_t digestBytes(std::uint64_t hash, const void* data,
                                 std::size_t n) {
  const unsigned char* bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i = 0; i < n; ++i) {
    hash ^= static_cast<std::uint64_t>(bytes[i]);
    hash *= kFnvPrime;
  }
  return hash;
}

template <typename U>
void writeValue(std::ostream& out, const U& value, std::uint64_t* digest) {
  if (digest != nullptr) *digest = digestBytes(*digest, &value, sizeof(U));
  out.write(reinterpret_cast<const char*>(&value), sizeof(U));
  if (!out) fail("write failed");
}

inline void writeBytes(std::ostream& out, const void* data, std::size_t n,
                       std::uint64_t* digest) {
  if (digest != nullptr) *digest = digestBytes(*digest, data, n);
  out.write(static_cast<const char*>(data), static_cast<std::streamsize>(n));
  if (!out) fail("write failed");
}

class Reader {
 public:
  Reader(std::istream& in, std::streamoff size) : _in(in), _remaining(size) {}

  std::streamoff remaining() const { return _remaining; }

  void readRaw(void* dst, std::size_t n, const char* what) {
    if (_remaining < 0 || static_cast<std::uint64_t>(_remaining) < n) {
      fail(std::string("truncated file while reading ") + what);
    }
    _in.read(static_cast<char*>(dst), static_cast<std::streamsize>(n));
    if (!_in) fail(std::string("read failed at ") + what);
    _remaining -= static_cast<std::streamoff>(n);
  }

  template <typename U>
  U value(const char* what) {
    U result{};
    readRaw(&result, sizeof(U), what);
    return result;
  }

  template <typename U>
  U hashed(std::uint64_t& digest, const char* what) {
    const U result = value<U>(what);
    digest = digestBytes(digest, &result, sizeof(U));
    return result;
  }

  void hashedBytes(void* dst, std::size_t n, std::uint64_t& digest,
                   const char* what) {
    readRaw(dst, n, what);
    digest = digestBytes(digest, dst, n);
  }

 private:
  std::istream& _in;
  std::streamoff _remaining;
};

inline std::uint64_t cellProduct(std::int32_t nx, std::int32_t ny,
                                 std::int32_t nz) {
  const std::uint64_t dims[3] = {static_cast<std::uint64_t>(nx),
                                 static_cast<std::uint64_t>(ny),
                                 static_cast<std::uint64_t>(nz)};
  std::uint64_t total = 1;
  for (std::uint64_t d : dims) {
    if (d == 0) fail("grid dimension must be positive");
    if (total > std::numeric_limits<std::uint64_t>::max() / d) {
      fail("grid size overflows");
    }
    total *= d;
  }
  return total;
}

inline void checkSortedCell(std::uint64_t cell, std::uint64_t n, bool hasPrevious,
                            std::uint64_t previous, const char* what) {
  if (cell >= n) fail(std::string(what) + " cell id is outside the grid");
  if (hasPrevious && cell <= previous) {
    fail(std::string(what) + " cells are not sorted and unique");
  }
}

}  // namespace detail

template <typename T>
void writeGeometryFile(const std::string& path, const VoxelGeometry<T>& geometry,
                       const GeometryRecipe<T>& recipe) {
  static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value,
                "geometry_file stores float or double only");

  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) detail::fail("cannot open '" + path + "' for writing");

  out.write(detail::kGeometryFileMagic, sizeof(detail::kGeometryFileMagic));
  if (!out) detail::fail("write failed");

  detail::writeValue(out, detail::kGeometryFileVersion, nullptr);
  detail::writeValue(out, std::is_same<T, float>::value ? std::uint32_t{0}
                                                        : std::uint32_t{1},
                     nullptr);
  detail::writeValue(out, std::uint32_t{3}, nullptr);
  detail::writeValue(out, static_cast<std::int32_t>(geometry.getNx()), nullptr);
  detail::writeValue(out, static_cast<std::int32_t>(geometry.getNy()), nullptr);
  detail::writeValue(out, static_cast<std::int32_t>(geometry.getNz()), nullptr);
  const Vector<T, 3>& min = geometry.getDomain().getMin();
  const Vector<T, 3>& max = geometry.getDomain().getMax();
  for (int i = 0; i < 3; ++i) detail::writeValue(out, min[i], nullptr);
  for (int i = 0; i < 3; ++i) detail::writeValue(out, max[i], nullptr);
  detail::writeValue(out, geometry.getVoxelSize(), nullptr);
  detail::writeValue(out, recipe.BandWidth, nullptr);

  detail::writeValue(out, static_cast<std::uint64_t>(geometry.size()), nullptr);
  std::uint64_t digest = detail::kFnvOffset;
  for (std::size_t id = 0; id < geometry.size(); ++id) {
    detail::writeValue(out, static_cast<std::uint8_t>(geometry.getRole(id)),
                       &digest);
  }
  for (std::size_t id = 0; id < geometry.size(); ++id) {
    detail::writeValue(out, geometry.getEntity(id), &digest);
  }
  for (std::size_t id = 0; id < geometry.size(); ++id) {
    detail::writeValue(out, geometry.getMaterial(id), &digest);
  }

  if (geometry.getPatchCount() > std::numeric_limits<std::uint32_t>::max()) {
    detail::fail("patch count does not fit the file format");
  }
  detail::writeValue(out, static_cast<std::uint32_t>(geometry.getPatchCount()),
                     &digest);
  for (std::size_t p = 0; p < geometry.getPatchCount(); ++p) {
    const BoundaryPatch& patch = geometry.getPatch(p);
    detail::writeValue(out, patch.Entity, &digest);
    detail::writeValue(out, static_cast<std::uint8_t>(patch.Kind), &digest);
    detail::writeValue(out, std::uint8_t{0}, &digest);
    detail::writeValue(out, static_cast<std::uint64_t>(patch.Cells.size()),
                       &digest);
    for (std::size_t id : patch.Cells) {
      detail::writeValue(out, static_cast<std::uint64_t>(id), &digest);
    }
  }

  detail::writeValue(out, static_cast<std::uint64_t>(geometry.getBandSize()),
                     &digest);
  for (std::size_t id : geometry.getBandCells()) {
    const SurfaceSample<T>& sample = geometry.getSurfaceSample(id);
    detail::writeValue(out, static_cast<std::uint64_t>(id), &digest);
    detail::writeValue(out, sample.Distance, &digest);
    for (int i = 0; i < 3; ++i) {
      detail::writeValue(out, sample.Normal[i], &digest);
    }
    detail::writeValue(out, sample.Entity, &digest);
    detail::writeValue(out, std::uint16_t{0}, &digest);
  }

  if (recipe.Sources.size() > std::numeric_limits<std::uint32_t>::max()) {
    detail::fail("source count does not fit the file format");
  }
  detail::writeValue(out, static_cast<std::uint32_t>(recipe.Sources.size()),
                     &digest);
  for (const GeometrySource<T>& source : recipe.Sources) {
    detail::writeValue(out, static_cast<std::uint8_t>(source.Kind), &digest);
    detail::writeValue(out, source.Entity, &digest);
    detail::writeValue(out, source.Material, &digest);
    detail::writeValue(out, static_cast<std::uint8_t>(source.Role), &digest);
    detail::writeValue(out, static_cast<std::uint8_t>(source.RequireClosed ? 1 : 0),
                       &digest);

    switch (source.Kind) {
      case GeometrySourceKind::Stl: {
        if (source.Path.size() > detail::kMaxPathLength) {
          detail::fail("source path exceeds " +
                       std::to_string(detail::kMaxPathLength) + " bytes");
        }
        detail::writeValue(out, static_cast<std::uint32_t>(source.Path.size()),
                           &digest);
        detail::writeBytes(out, source.Path.data(), source.Path.size(), &digest);
        break;
      }
      case GeometrySourceKind::Box: {
        const Vector<T, 3>& boxMin = source.Box.getMin();
        const Vector<T, 3>& boxMax = source.Box.getMax();
        for (int i = 0; i < 3; ++i) detail::writeValue(out, boxMin[i], &digest);
        for (int i = 0; i < 3; ++i) detail::writeValue(out, boxMax[i], &digest);
        break;
      }
      case GeometrySourceKind::Sphere: {
        for (int i = 0; i < 3; ++i) {
          detail::writeValue(out, source.Center[i], &digest);
        }
        detail::writeValue(out, source.Radius, &digest);
        break;
      }
    }
  }

  detail::writeValue(out, digest, nullptr);
}

template <typename T>
GeometryFile<T> GeometryFile<T>::readFromFile(const std::string& path) {
  static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value,
                "geometry_file stores float or double only");

  std::ifstream in(path, std::ios::binary);
  if (!in) detail::fail("cannot open '" + path + "' for reading");
  in.seekg(0, std::ios::end);
  const std::streamoff fileSize = in.tellg();
  if (fileSize < 0) detail::fail("cannot determine the size of '" + path + "'");
  in.seekg(0, std::ios::beg);
  detail::Reader reader(in, fileSize);

  char magic[sizeof(detail::kGeometryFileMagic)];
  reader.readRaw(magic, sizeof(magic), "magic");
  if (std::memcmp(magic, detail::kGeometryFileMagic, sizeof(magic)) != 0) {
    detail::fail("bad magic, not a FLBGEO file");
  }

  const std::uint32_t version = reader.value<std::uint32_t>("version");
  if (version != detail::kGeometryFileVersion) {
    detail::fail("unsupported version " + std::to_string(version));
  }
  const std::uint32_t scalarTag = reader.value<std::uint32_t>("scalar tag");
  const std::uint32_t expectedScalarTag =
    std::is_same<T, float>::value ? 0u : 1u;
  if (scalarTag != expectedScalarTag) {
    detail::fail("scalar type tag does not match this build");
  }
  const std::uint32_t dim = reader.value<std::uint32_t>("dimension");
  if (dim != 3) detail::fail("only 3D geometry files are supported");

  const std::int32_t nx = reader.value<std::int32_t>("nx");
  const std::int32_t ny = reader.value<std::int32_t>("ny");
  const std::int32_t nz = reader.value<std::int32_t>("nz");
  if (nx <= 0 || ny <= 0 || nz <= 0) {
    detail::fail("grid dimensions must be positive");
  }

  T min[3];
  T max[3];
  for (int i = 0; i < 3; ++i) min[i] = reader.value<T>("domain min");
  for (int i = 0; i < 3; ++i) max[i] = reader.value<T>("domain max");
  const T voxelSize = reader.value<T>("voxel size");
  const T bandWidth = reader.value<T>("band width");

  const std::uint64_t cellCount = reader.value<std::uint64_t>("cell count");
  if (cellCount != detail::cellProduct(nx, ny, nz)) {
    detail::fail("cell count does not match the grid");
  }
  if (cellCount > std::numeric_limits<std::size_t>::max()) {
    detail::fail("cell count does not fit in memory");
  }
  const std::size_t cells = static_cast<std::size_t>(cellCount);
  if (cells > std::numeric_limits<std::size_t>::max() / sizeof(std::uint16_t)) {
    detail::fail("cell count is too large");
  }
  // A matching digest cannot be checked before the body is read, but at least
  // reject a header that claims more cells than the file could hold.
  constexpr std::uint64_t kFixedTail =
    4 + 8 + 4 + 8;  // patchCount, bandCount, sourceCount, digest
  if (cellCount >
      (std::numeric_limits<std::uint64_t>::max() - kFixedTail) / 5) {
    detail::fail("cell count overflows the file size check");
  }
  const std::uint64_t minimumBody = cellCount * 5 + kFixedTail;
  if (minimumBody > static_cast<std::uint64_t>(reader.remaining())) {
    detail::fail("file is too small for the declared grid");
  }

  const Vector<T, 3> domainMin{min[0], min[1], min[2]};
  const Vector<T, 3> domainMax{max[0], max[1], max[2]};
  GeometryRecipe<T> recipe;
  recipe.BandWidth = bandWidth;
  GeometryFile<T> file{
    VoxelGeometry<T>(nx, ny, nz, AABB<T, 3>(domainMin, domainMax), voxelSize),
    std::move(recipe)};
  VoxelGeometry<T>& geometry = file.Geometry;

  std::uint64_t digest = detail::kFnvOffset;
  std::vector<std::uint8_t> roles(cells);
  reader.hashedBytes(roles.data(), roles.size(), digest, "role array");
  for (std::uint8_t role : roles) {
    if (role > static_cast<std::uint8_t>(VoxelRole::Solid)) {
      detail::fail("invalid role value in the role array");
    }
  }

  reader.hashedBytes(geometry._entity.data(),
                     cells * sizeof(std::uint16_t), digest, "entity array");
  reader.hashedBytes(geometry._material.data(),
                     cells * sizeof(std::uint16_t), digest, "material array");
  for (std::size_t id = 0; id < cells; ++id) {
    geometry._role[id] = static_cast<VoxelRole>(roles[id]);
  }

  const std::uint32_t patchCount =
    reader.hashed<std::uint32_t>(digest, "patch count");
  if (patchCount > static_cast<std::uint64_t>(reader.remaining()) / 12) {
    detail::fail("patch count exceeds the remaining file size");
  }
  geometry._patches.resize(patchCount);
  for (std::uint32_t p = 0; p < patchCount; ++p) {
    BoundaryPatch patch;
    patch.Entity = reader.hashed<std::uint16_t>(digest, "patch entity");
    const std::uint8_t kind = reader.hashed<std::uint8_t>(digest, "patch kind");
    const std::uint8_t reserved =
      reader.hashed<std::uint8_t>(digest, "patch reserved");
    (void)reserved;
    if (kind > static_cast<std::uint8_t>(PatchKind::Symmetry)) {
      detail::fail("invalid patch kind");
    }
    patch.Kind = static_cast<PatchKind>(kind);

    const std::uint64_t count =
      reader.hashed<std::uint64_t>(digest, "patch cell count");
    if (count > cellCount) detail::fail("patch cell count exceeds the grid");
    patch.Cells.resize(static_cast<std::size_t>(count));
    std::uint64_t previous = 0;
    for (std::size_t i = 0; i < patch.Cells.size(); ++i) {
      const std::uint64_t cell =
        reader.hashed<std::uint64_t>(digest, "patch cell");
      detail::checkSortedCell(cell, cellCount, i != 0, previous, "patch");
      previous = cell;
      patch.Cells[i] = static_cast<std::size_t>(cell);
    }
    geometry._patches[p] = std::move(patch);
  }

  const std::uint64_t bandCount =
    reader.hashed<std::uint64_t>(digest, "band count");
  if (bandCount > cellCount) detail::fail("band count exceeds the grid");
  geometry._bandIds.resize(static_cast<std::size_t>(bandCount));
  geometry._band.reserve(static_cast<std::size_t>(bandCount));
  std::uint64_t previousBand = 0;
  for (std::size_t i = 0; i < geometry._bandIds.size(); ++i) {
    const std::uint64_t cell =
      reader.hashed<std::uint64_t>(digest, "band cell");
    detail::checkSortedCell(cell, cellCount, i != 0, previousBand, "band");
    previousBand = cell;
    const T distance = reader.hashed<T>(digest, "band distance");
    T normal[3];
    for (int j = 0; j < 3; ++j) normal[j] = reader.hashed<T>(digest, "band normal");
    const std::uint16_t entity =
      reader.hashed<std::uint16_t>(digest, "band entity");
    const std::uint16_t reserved =
      reader.hashed<std::uint16_t>(digest, "band reserved");
    (void)reserved;

    geometry._bandIds[i] = static_cast<std::size_t>(cell);
    geometry._band.emplace(
      static_cast<std::size_t>(cell),
      SurfaceSample<T>{distance, Vector<T, 3>{normal[0], normal[1], normal[2]},
                       entity});
  }

  const std::uint32_t sourceCount =
    reader.hashed<std::uint32_t>(digest, "source count");
  // Fixed source prefix (7 bytes) plus the smallest body (an empty STL path,
  // 4 bytes) is 11 bytes.
  if (sourceCount > static_cast<std::uint64_t>(reader.remaining()) / 11) {
    detail::fail("source count exceeds the remaining file size");
  }
  file.Recipe.Sources.resize(sourceCount);
  for (std::uint32_t s = 0; s < sourceCount; ++s) {
    const std::uint8_t kind =
      reader.hashed<std::uint8_t>(digest, "source kind");
    const std::uint16_t entity =
      reader.hashed<std::uint16_t>(digest, "source entity");
    const std::uint16_t material =
      reader.hashed<std::uint16_t>(digest, "source material");
    const std::uint8_t role = reader.hashed<std::uint8_t>(digest, "source role");
    const std::uint8_t flags = reader.hashed<std::uint8_t>(digest, "source flags");
    if (kind > static_cast<std::uint8_t>(GeometrySourceKind::Sphere)) {
      detail::fail("invalid source kind");
    }
    if (role > static_cast<std::uint8_t>(VoxelRole::Solid)) {
      detail::fail("invalid source role");
    }

    GeometrySource<T> source;
    source.Kind = static_cast<GeometrySourceKind>(kind);
    source.Entity = entity;
    source.Material = material;
    source.Role = static_cast<VoxelRole>(role);
    source.RequireClosed = (flags & std::uint8_t{1}) != 0;

    switch (source.Kind) {
      case GeometrySourceKind::Stl: {
        const std::uint32_t length =
          reader.hashed<std::uint32_t>(digest, "source path length");
        if (length > detail::kMaxPathLength) {
          detail::fail("source path is longer than " +
                       std::to_string(detail::kMaxPathLength) + " bytes");
        }
        source.Path.resize(length);
        reader.hashedBytes(source.Path.data(), source.Path.size(), digest,
                           "source path");
        break;
      }
      case GeometrySourceKind::Box: {
        T boxMin[3];
        T boxMax[3];
        for (int i = 0; i < 3; ++i) boxMin[i] = reader.hashed<T>(digest, "box min");
        for (int i = 0; i < 3; ++i) boxMax[i] = reader.hashed<T>(digest, "box max");
        source.Box = AABB<T, 3>(Vector<T, 3>{boxMin[0], boxMin[1], boxMin[2]},
                                Vector<T, 3>{boxMax[0], boxMax[1], boxMax[2]});
        break;
      }
      case GeometrySourceKind::Sphere: {
        T center[3];
        for (int i = 0; i < 3; ++i) {
          center[i] = reader.hashed<T>(digest, "sphere center");
        }
        source.Center = Vector<T, 3>{center[0], center[1], center[2]};
        source.Radius = reader.hashed<T>(digest, "sphere radius");
        break;
      }
    }
    file.Recipe.Sources[s] = std::move(source);
  }

  const std::uint64_t storedDigest =
    reader.value<std::uint64_t>("geometry digest");
  if (storedDigest != digest) detail::fail("geometry digest mismatch");
  if (reader.remaining() != 0) detail::fail("trailing bytes after the digest");

  return file;
}

template <typename T>
GeometryFile<T> readGeometryFile(const std::string& path) {
  return GeometryFile<T>::readFromFile(path);
}

}  // namespace flbgeo
