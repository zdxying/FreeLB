/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// checkpoint.h
// declarations for checkpoint / restart
//
// Key points:
//   - only "fields + step count" are serialized; the geometry is replayed by
//     the example's own setup code (same parameters assumed);
//   - POP is exported through the container's **logical view**; the file holds
//     no rotation state, which decouples it from the POP container;
//   - every field manager (block-lattice or plain block-field) is written as a
//     whole, in its field-pack order;
//   - checkCompatible() validates the block table per block before loading;
//     any mismatch is a hard error, never silent;
//   - under CUDA the module does copyToHost before writing and copyToDevice
//     after reading.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// BlockFieldManager / BlockField
#include "data_struct/field_struct.h"
// BlockLatticeManagerBase: a block-lattice manager is accepted directly; its
// whole field pack (owned + referenced) is written together
#include "data_struct/block_lattice_base.h"
// MPI_RANK / mpi().barrier()
#include "parallel/mpi_manager.h"
// DirCreator: creates the snapshot directory
#include "utils/directories.h"
// Vector<T,D>: elements may be vectors, written as D*sizeof(scalar) raw bytes
#include "data_struct/Vector.h"

namespace flbckpt {

// element type tag stored in the file.  Vector<T,D> records its scalar type;
// the byte count is D * sizeof(scalar).
enum class ElemTag : std::uint32_t { f32 = 0, f64 = 1, u8 = 2, u16 = 3, half = 4 };

// format version; bump on an incompatible change, old files are rejected
constexpr std::uint32_t kFormatVersion = 1;
// binary file magic
constexpr const char* kMagic = "FLBCKPT";

// geometric identity of one block, validated block by block on load
struct BlockRecord {
  std::int32_t blockid = -1;
  std::uint8_t level = 0;
  std::int32_t overlap = 0;
  std::uint32_t dim = 0;
  std::vector<std::int32_t> mesh;
  std::uint64_t n = 0;
};

// metadata of one field (counts[b][j] = length of array j of block b)
struct FieldRecord {
  std::uint32_t index = 0;
  std::uint32_t array_dim = 0;
  ElemTag tag = ElemTag::f32;
  std::uint32_t elem_dim = 1;  // 1 = scalar, 2/3 = Vector
  std::string name;
  std::vector<std::vector<std::uint64_t>> counts;
};

// file header + all metadata (no payload)
struct CheckpointHeader {
  std::uint32_t version = 0;
  std::int64_t step = 0;
  std::int32_t rank = 0;
  std::int32_t nranks = 1;
  std::uint32_t dim = 0;
  std::vector<BlockRecord> blocks;
  std::vector<FieldRecord> fields;
};

// ---------------------------------------------------------------------------
// public API: variadic, accepts any number of "block-lattice managers"
// (BlockLatticeManager and its derivatives, CA included)
// and "block-field managers" (e.g. FlagFM).  Fields are written whole, in
// their own field-pack order.
// ---------------------------------------------------------------------------

// write a snapshot at a quiescent point (communication done, before the next Stream)
template <typename... FMs>
void write(const std::string& dir, std::int64_t step, FMs&... fms);

// restore; step < 0 reads the "latest" pointer file.  Returns the restored step
template <typename... FMs>
std::int64_t read(const std::string& dir, std::int64_t step, FMs&... fms);

// validate without loading: block table / element type / array sizes against
// the current geometry
template <typename... FMs>
bool checkCompatible(const std::string& dir, std::int64_t step, FMs&... fms);

// read a snapshot header (no payload); returns false on failure
bool readHeader(const std::string& dir, std::int64_t step, CheckpointHeader& hdr);

// step of the newest complete snapshot in the directory, or -1 if there is none
std::int64_t latestStep(const std::string& dir);

}  // namespace flbckpt
