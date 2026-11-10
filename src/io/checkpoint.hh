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

// checkpoint.hh
// implementation of checkpoint / restart
//
// file layout (one binary per rank):
//   magic "FLBCKPT" | u32 version | i64 step | i32 rank | i32 nranks | u32 dim
//   u64 nblocks { i32 blockid | u8 level | i32 overlap | u32 dim | u32 nmesh |
//                 i32 mesh[nmesh] | u64 n }
//   u64 nfields { u32 index | u32 array_dim | u32 tag | u32 elem_dim | u32 namelen |
//                 char name[] | u64 nblocks { u32 array_dim | u64 counts[] } }
//   payload: in meta field order, per field per block per array, u64 count +
//            raw bytes
//        (a Data<> scalar field has no count prefix; its size comes from
//         tag/elem_dim, and meta marks such fields with array_dim == 0)
//
// Restore assumption: replaying the setup with the same parameters gives a
// block-for-block identical geometry, so only fields and the step count are
// restored; a block-table mismatch is a hard error.

#include "io/checkpoint.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
// POSIX directory listing (used to discover rank*.ck files)
#include <dirent.h>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#if defined(FREELB_POP_FP16) && defined(__CUDACC__)
#include <cuda_fp16.h>
#endif

namespace flbckpt {
namespace detail {

// ---------------------------------------------------------------------------
// low-level I/O
//
// Note: always use IF_MPI_RANK, never MPI_RANK -- the latter expands to
// if (rank != 0) { return; }, which returns early from a loop or a template
// function.
// ---------------------------------------------------------------------------
template <typename U>
void wVal(std::ostream& os, const U& v) {
  os.write(reinterpret_cast<const char*>(&v), static_cast<std::streamsize>(sizeof(U)));
}
template <typename U>
void rVal(std::istream& is, U& v) {
  is.read(reinterpret_cast<char*>(&v), static_cast<std::streamsize>(sizeof(U)));
}
inline void wStr(std::ostream& os, const std::string& s) {
  wVal<std::uint32_t>(os, static_cast<std::uint32_t>(s.size()));
  if (!s.empty()) os.write(s.data(), static_cast<std::streamsize>(s.size()));
}
inline std::string rStr(std::istream& is) {
  std::uint32_t n = 0;
  rVal(is, n);
  std::string s(static_cast<std::size_t>(n), '\0');
  if (n) is.read(&s[0], static_cast<std::streamsize>(n));
  return s;
}

[[noreturn]] inline void die(const std::string& msg) {
  IF_MPI_RANK(0) { std::cerr << "[flbckpt] " << msg << std::endl; }
  // exit on all ranks together, so the others do not block in communication
  exit(1);
}

// ---------------------------------------------------------------------------
// element type: Vector<T,D> records its scalar type and dimension, and is
// written as D*sizeof(scalar) bytes
// ---------------------------------------------------------------------------
template <typename V>
struct ElemTraits {
  using scalar = V;
  static constexpr unsigned int dim = 1;
};
template <typename T, unsigned int D>
struct ElemTraits<Vector<T, D>> {
  using scalar = T;
  static constexpr unsigned int dim = D;
};

template <typename V>
constexpr ElemTag tagOfScalar() {
  if constexpr (std::is_same_v<V, float>) {
    return ElemTag::f32;
  } else if constexpr (std::is_same_v<V, double>) {
    return ElemTag::f64;
  } else if constexpr (std::is_same_v<V, std::uint8_t>) {
    return ElemTag::u8;
  } else if constexpr (std::is_same_v<V, std::uint16_t>) {
    return ElemTag::u16;
#if defined(FREELB_POP_FP16) && defined(__CUDACC__)
  } else if constexpr (std::is_same_v<V, __half>) {
    return ElemTag::half;
#endif
  } else {
    static_assert(sizeof(V) == 0, "unsupported field element type for a checkpoint");
  }
}

template <typename FieldType>
ElemTag fieldTag() {
  return tagOfScalar<typename ElemTraits<typename FieldType::value_type>::scalar>();
}
template <typename FieldType>
constexpr unsigned int fieldElemDim() {
  return ElemTraits<typename FieldType::value_type>::dim;
}

inline std::size_t elemSizeOf(const FieldRecord& r) {
  const std::size_t s = r.tag == ElemTag::f32   ? 4
                        : r.tag == ElemTag::f64 ? 8
                        : r.tag == ElemTag::u8  ? 1
                        : r.tag == ElemTag::u16 ? 2
                                                : 2 /* half */;
  return s * r.elem_dim;
}

inline const char* tagName(ElemTag t) {
  switch (t) {
    case ElemTag::f32: return "f32";
    case ElemTag::f64: return "f64";
    case ElemTag::u8: return "u8";
    case ElemTag::u16: return "u16";
    case ElemTag::half: return "half";
  }
  return "?";
}

inline std::string meshToString(const std::vector<std::int32_t>& mesh) {
  std::ostringstream os;
  for (std::size_t i = 0; i < mesh.size(); ++i) os << (i ? "," : "{") << mesh[i];
  os << "}";
  return os.str();
}

// ---------------------------------------------------------------------------
// logical-view read/write
//
// The POP container (CyclicArray) stores populations rotated, so its logical
// view is two disjoint ranges in memory; write and read range by range
// instead of dumping the raw buffer -- the file holds no rotation state.
// ---------------------------------------------------------------------------
template <typename Value, typename ArrayType>
void writeLogical(std::ostream& os, const ArrayType& arr, std::size_t n) {
  if (n == 0) return;
  const Value* p = arr.getdataPtr(0);
  // fast path: the logical view is contiguous (GenericArray, and a CyclicArray
  // whose shift has been normalized)
  if (arr.getdataPtr(n - 1) == p + static_cast<std::ptrdiff_t>(n - 1)) {
    os.write(reinterpret_cast<const char*>(p),
             static_cast<std::streamsize>(n * sizeof(Value)));
    return;
  }
  std::size_t i = 0;
  while (i < n) {
    p = arr.getdataPtr(i);
    std::size_t j = i + 1;
    while (j < n && arr.getdataPtr(j) == p + static_cast<std::ptrdiff_t>(j - i)) ++j;
    os.write(reinterpret_cast<const char*>(p),
             static_cast<std::streamsize>((j - i) * sizeof(Value)));
    i = j;
  }
}

template <typename Value, typename ArrayType>
void readLogical(std::istream& is, ArrayType& arr, std::size_t n) {
  if (n == 0) return;
  Value* p = arr.getdataPtr(0);
  if (arr.getdataPtr(n - 1) == p + static_cast<std::ptrdiff_t>(n - 1)) {
    is.read(reinterpret_cast<char*>(p), static_cast<std::streamsize>(n * sizeof(Value)));
    return;
  }
  std::size_t i = 0;
  while (i < n) {
    p = arr.getdataPtr(i);
    std::size_t j = i + 1;
    while (j < n && arr.getdataPtr(j) == p + static_cast<std::ptrdiff_t>(j - i)) ++j;
    is.read(reinterpret_cast<char*>(p), static_cast<std::streamsize>((j - i) * sizeof(Value)));
    i = j;
  }
}

// does the container provide normalize() (reset the rotation to the unrotated
// view)?  CyclicArray / StreamMapArray do
template <typename A, typename = void>
struct HasNormalize : std::false_type {};
template <typename A>
struct HasNormalize<A, std::void_t<decltype(std::declval<A&>().normalize())>>
    : std::true_type {};

// does the container provide refreshDeviceObject() (re-upload the device
// mirror object, see field.h)?
template <typename A, typename = void>
struct HasRefreshDeviceObject : std::false_type {};
template <typename A>
struct HasRefreshDeviceObject<
  A, std::void_t<decltype(std::declval<A&>().refreshDeviceObject())>> : std::true_type {};

// ---------------------------------------------------------------------------
// write/restore one field (the whole manager: every block)
// ---------------------------------------------------------------------------
template <typename FieldType, typename T, unsigned int Dim>
void dumpField(std::ostream& os, const BlockFieldManager<FieldType, T, Dim>& fm) {
  using Value = typename FieldType::value_type;
  for (const auto& bf : fm.getBlockFields()) {
    if constexpr (FieldType::isField) {
      for (unsigned int j = 0; j < FieldType::array_dim; ++j) {
        const auto& arr = bf.getField(j);
        const std::size_t n = arr.size();
        wVal<std::uint64_t>(os, static_cast<std::uint64_t>(n));
        writeLogical<Value>(os, arr, n);
      }
    } else {
      // Data<T,Base>: one scalar per block, no count prefix
      const Value& v = bf.get();
      os.write(reinterpret_cast<const char*>(&v),
               static_cast<std::streamsize>(sizeof(Value)));
    }
  }
}

template <typename FieldType, typename T, unsigned int Dim>
// payload size of **one block** in one field (used to skip blocks owned by
// another rank)
inline std::uint64_t blockPayloadSize(const FieldRecord& rec, std::size_t b) {
  if (rec.array_dim == 0) return static_cast<std::uint64_t>(elemSizeOf(rec));
  std::uint64_t sz = 0;
  for (std::uint64_t v : rec.counts[b])
    sz += static_cast<std::uint64_t>(sizeof(std::uint64_t)) + v * elemSizeOf(rec);
  return sz;
}

// read the payload of **one block** (all arrays of the field, or the Data<>
// scalar)
template <typename FieldType, typename T, unsigned int Dim>
void loadOneBlock(std::istream& is, BlockField<FieldType, T, Dim>& bf,
                  const FieldRecord& rec) {
  using Value = typename FieldType::value_type;
  if constexpr (FieldType::isField) {
    for (unsigned int j = 0; j < FieldType::array_dim; ++j) {
      auto& arr = bf.getField(j);
      std::uint64_t n = 0;
      rVal(is, n);
      if (n != static_cast<std::uint64_t>(arr.size())) {
        std::ostringstream os;
        os << "field '" << rec.name << "' array " << j << ": file has " << n
           << " elements, the live lattice has " << arr.size()
           << " (incompatible geometry?)";
        die(os.str());
      }
      // Reset the rotation to the unrotated view first, then read through the
      // logical view.  The order matters:
      // reading first and normalizing afterwards writes the data at the logical
      // positions of the old rotation, and normalize() then redefines the
      // logical view as the physical order -- the data gets permuted.
      // After normalization the evolution depends only on the logical view and
      // Offset (Offset is reset by the example's Init() and does not change
      // inside the loop), so it is equivalent to an uninterrupted run.
      if constexpr (HasNormalize<std::decay_t<decltype(arr)>>::value) {
        arr.normalize();
#ifdef __CUDACC__
        // The device object caches start[] inside the object and copyToDevice()
        // does not push it, so it has to be re-uploaded after normalize(),
        // otherwise the device side reads a shifted window.
        if constexpr (HasRefreshDeviceObject<std::decay_t<decltype(arr)>>::value) {
          arr.refreshDeviceObject();
        }
#endif
      }
      readLogical<Value>(is, arr, static_cast<std::size_t>(n));
    }
  } else {
    Value& v = bf.get();
    is.read(reinterpret_cast<char*>(&v), static_cast<std::streamsize>(sizeof(Value)));
  }
}

// skip the payload of **one block** in one field
inline void skipOneBlock(std::istream& is, const FieldRecord& rec, std::size_t b) {
  if (rec.array_dim == 0) {
    // Data<> scalar field: one value per block, no count prefix
    is.seekg(static_cast<std::streamoff>(elemSizeOf(rec)), std::ios::cur);
    return;
  }
  for (std::uint64_t v : rec.counts[b])
    is.seekg(
      static_cast<std::streamoff>(sizeof(std::uint64_t) + v * elemSizeOf(rec)),
      std::ios::cur);
}

// skip the payload of every block of one field (used when this build does not
// have the field)
inline void skipField(std::istream& is, const FieldRecord& rec) {
  for (std::size_t b = 0; b < rec.counts.size(); ++b) skipOneBlock(is, rec, b);
}

// ---------------------------------------------------------------------------
// metadata collection
// ---------------------------------------------------------------------------
template <typename FieldType, typename T, unsigned int Dim>
void blockRecords(const BlockFieldManager<FieldType, T, Dim>& fm,
                  std::vector<BlockRecord>& recs) {
  for (const auto& bf : fm.getBlockFields()) {
    const auto& blk = bf.getBlock();
    BlockRecord r;
    r.blockid = blk.getBlockId();
    r.level = blk.getLevel();
    r.overlap = blk.getOverlap();
    r.dim = Dim;
    const Vector<int, Dim>& mesh = blk.getMesh();
    for (unsigned int i = 0; i < Dim; ++i) r.mesh.push_back(mesh[i]);
    r.n = static_cast<std::uint64_t>(blk.getN());
    recs.push_back(std::move(r));
  }
}

template <typename FieldType, typename T, unsigned int Dim>
void fieldCounts(const BlockFieldManager<FieldType, T, Dim>& fm,
                 std::vector<std::vector<std::uint64_t>>& counts) {
  counts.clear();
  if constexpr (FieldType::isField) {
    for (const auto& bf : fm.getBlockFields()) {
      std::vector<std::uint64_t> c(FieldType::array_dim);
      for (unsigned int j = 0; j < FieldType::array_dim; ++j)
        c[j] = static_cast<std::uint64_t>(bf.getField(j).size());
      counts.push_back(std::move(c));
    }
  } else {
    counts.assign(fm.size(), std::vector<std::uint64_t>{});
  }
}

// detect a "block-lattice manager" (it has FIELDS/FIELDPTRS member types);
// otherwise treat it as a plain block-field manager
template <typename MGR, typename = void>
struct IsBlockLatManager : std::false_type {};
template <typename MGR>
struct IsBlockLatManager<MGR, std::void_t<typename MGR::FIELDS, typename MGR::FIELDPTRS>>
    : std::true_type {};

// visit every field in a field pack: func(block field manager)
template <typename MGR, typename... Fields, typename Func>
void forEachFieldInPack(MGR& mgr, TypePack<Fields...>, Func&& func) {
  (func(mgr.template getField<Fields>()), ...);
}

// collect the live metadata and register the payload write/read function of
// every field
struct MetaCollector {
  enum class Mode { Write, Read, Check };

  CheckpointHeader live;
  Mode mode = Mode::Write;
  std::uint32_t fieldIndex = 0;
  bool haveBlocks = false;
  // the write side still writes a whole manager in order; the read side works
  // block by block, because a cross-rank restore picks only its own blocks
  // out of every rank's file (the others are skipped with seek)
  using PayloadWriter = std::function<void(std::ostream&)>;
  using BlockReader = std::function<void(std::istream&, std::size_t, const FieldRecord&)>;
  std::vector<std::pair<std::uint32_t, PayloadWriter>> writers;
  std::vector<std::pair<std::uint32_t, BlockReader>> readers;
  // read path: push the data back to the device after the payload is loaded
  std::vector<std::function<void()>> postRead;

  template <typename FM>
  void run(FM& fm) {
    if constexpr (IsBlockLatManager<std::decay_t<FM>>::value) {
      forEachFieldInPack(fm, typename std::decay_t<FM>::FIELDS{}, *this);
      forEachFieldInPack(fm, typename std::decay_t<FM>::FIELDPTRS{}, *this);
    } else {
      (*this)(fm);
    }
  }

  template <typename FieldType, typename T, unsigned int Dim>
  void operator()(BlockFieldManager<FieldType, T, Dim>& fm) {
    if (!haveBlocks) {
      blockRecords(fm, live.blocks);
      live.dim = Dim;
      haveBlocks = true;
    } else {
      // every manager shares the same geometry: the block tables must match
      std::vector<BlockRecord> recs;
      blockRecords(fm, recs);
      if (recs.size() != live.blocks.size())
        die("field managers disagree on the block count (" +
            std::to_string(recs.size()) + " vs " + std::to_string(live.blocks.size()) +
            "); only managers built from the same geometry can be checkpointed together");
      for (std::size_t i = 0; i < recs.size(); ++i) {
        if (recs[i].blockid != live.blocks[i].blockid ||
            recs[i].level != live.blocks[i].level || recs[i].n != live.blocks[i].n)
          die("field managers disagree on block " + std::to_string(i) +
              "; only managers built from the same geometry can be checkpointed together");
      }
    }

    FieldRecord rec;
    rec.index = fieldIndex;
    rec.array_dim = FieldType::isField ? FieldType::array_dim : 0;
    rec.tag = fieldTag<FieldType>();
    rec.elem_dim = fieldElemDim<FieldType>();
    rec.name = typeid(FieldType).name();
    fieldCounts(fm, rec.counts);
    live.fields.push_back(rec);

#ifdef __CUDACC__
    // write: copy the device data back to the host first; read: push it to the
    // device after the payload is loaded (postRead)
    if (mode == Mode::Write) fm.copyToHost();
    if (mode == Mode::Read) postRead.push_back([&fm]() { fm.copyToDevice(); });
#endif

    writers.emplace_back(fieldIndex, [&fm](std::ostream& os) { dumpField(os, fm); });
    readers.emplace_back(fieldIndex,
                         [&fm](std::istream& is, std::size_t b, const FieldRecord& r) {
                           loadOneBlock(is, fm.getBlockFields()[b], r);
                         });
    ++fieldIndex;
  }
};

}  // namespace detail

// ---------------------------------------------------------------------------
// file header
// ---------------------------------------------------------------------------
namespace detail {

void writeHeader(std::ostream& os, const CheckpointHeader& h) {
  os.write(kMagic, static_cast<std::streamsize>(std::strlen(kMagic)));
  wVal<std::uint32_t>(os, kFormatVersion);
  wVal<std::int64_t>(os, h.step);
  wVal<std::int32_t>(os, h.rank);
  wVal<std::int32_t>(os, h.nranks);
  wVal<std::uint32_t>(os, h.dim);
  wVal<std::uint64_t>(os, static_cast<std::uint64_t>(h.blocks.size()));
  for (const auto& b : h.blocks) {
    wVal<std::int32_t>(os, b.blockid);
    wVal<std::uint8_t>(os, b.level);
    wVal<std::int32_t>(os, b.overlap);
    wVal<std::uint32_t>(os, b.dim);
    wVal<std::uint32_t>(os, static_cast<std::uint32_t>(b.mesh.size()));
    for (int m : b.mesh) wVal<std::int32_t>(os, m);
    wVal<std::uint64_t>(os, b.n);
  }
  wVal<std::uint64_t>(os, static_cast<std::uint64_t>(h.fields.size()));
  for (const auto& f : h.fields) {
    wVal<std::uint32_t>(os, f.index);
    wVal<std::uint32_t>(os, f.array_dim);
    wVal<std::uint32_t>(os, static_cast<std::uint32_t>(f.tag));
    wVal<std::uint32_t>(os, f.elem_dim);
    wStr(os, f.name);
    wVal<std::uint64_t>(os, static_cast<std::uint64_t>(f.counts.size()));
    for (const auto& c : f.counts) {
      wVal<std::uint32_t>(os, static_cast<std::uint32_t>(c.size()));
      for (std::uint64_t v : c) wVal<std::uint64_t>(os, v);
    }
  }
}

void readHeaderFrom(std::istream& is, CheckpointHeader& h) {
  char magic[8] = {0};
  is.read(magic, static_cast<std::streamsize>(std::strlen(kMagic)));
  if (std::memcmp(magic, kMagic, std::strlen(kMagic)) != 0)
    die("not a FreeLB checkpoint file (bad magic)");
  rVal<std::uint32_t>(is, h.version);
  rVal<std::int64_t>(is, h.step);
  rVal<std::int32_t>(is, h.rank);
  rVal<std::int32_t>(is, h.nranks);
  rVal<std::uint32_t>(is, h.dim);
  if (h.version != kFormatVersion)
    die("checkpoint format version " + std::to_string(h.version) +
        " is not supported by this build (expects " + std::to_string(kFormatVersion) +
        "); the file must be re-created");
  std::uint64_t nb = 0;
  rVal<std::uint64_t>(is, nb);
  h.blocks.resize(static_cast<std::size_t>(nb));
  for (auto& b : h.blocks) {
    rVal<std::int32_t>(is, b.blockid);
    rVal<std::uint8_t>(is, b.level);
    rVal<std::int32_t>(is, b.overlap);
    rVal<std::uint32_t>(is, b.dim);
    std::uint32_t nm = 0;
    rVal<std::uint32_t>(is, nm);
    b.mesh.resize(nm);
    for (auto& m : b.mesh) rVal<std::int32_t>(is, m);
    rVal<std::uint64_t>(is, b.n);
  }
  std::uint64_t nf = 0;
  rVal<std::uint64_t>(is, nf);
  h.fields.resize(static_cast<std::size_t>(nf));
  for (auto& f : h.fields) {
    rVal<std::uint32_t>(is, f.index);
    rVal<std::uint32_t>(is, f.array_dim);
    std::uint32_t tag = 0;
    rVal<std::uint32_t>(is, tag);
    f.tag = static_cast<ElemTag>(tag);
    rVal<std::uint32_t>(is, f.elem_dim);
    f.name = rStr(is);
    std::uint64_t nblk = 0;
    rVal<std::uint64_t>(is, nblk);
    f.counts.resize(static_cast<std::size_t>(nblk));
    for (auto& c : f.counts) {
      std::uint32_t na = 0;
      rVal<std::uint32_t>(is, na);
      c.resize(na);
      for (auto& v : c) rVal<std::uint64_t>(is, v);
    }
  }
}

// ---------------------------------------------------------------------------
// cross-rank restore: find every rank file of this step and build a mapping by
// blockid
//
// blockid is a global id assigned when the blocks are created; LoadBalancing
// only changes ownership, not the id, so it can be used to match blocks
// across different decompositions.
// ---------------------------------------------------------------------------
// NB: no user-declared destructor -- it would suppress the implicit move
// constructor, and std::ifstream is not copyable, so returning this by value
// would then try to copy.  The streams close themselves on destruction.
struct RankFileSet {
  std::vector<std::string> paths;
  std::vector<std::ifstream> streams;
  std::vector<CheckpointHeader> headers;
  std::int64_t step = 0;
  std::int32_t nranks = 0;
};

inline RankFileSet openRankFiles(const std::string& dir, std::int64_t step,
                                 bool needPayload = true) {
  RankFileSet set;
  set.step = step;
  const std::string stepdir = dir + "/step_" + std::to_string(step);
  // list rank*.ck under stepdir
  std::vector<std::string> paths;
  if (DIR* d = opendir(stepdir.c_str())) {
    while (dirent* e = readdir(d)) {
      const std::string name = e->d_name;
      if (name.size() > 3 && name.compare(0, 4, "rank") == 0 &&
          name.compare(name.size() - 3, 3, ".ck") == 0)
        paths.push_back(stepdir + "/" + name);
    }
    closedir(d);
  }
  if (paths.empty()) die("no rank*.ck checkpoint files in " + stepdir);
  std::sort(paths.begin(), paths.end());
  const std::size_t nf = paths.size();
  set.paths = paths;
  set.streams.reserve(nf);
  set.headers.resize(nf);
  set.nranks = static_cast<std::int32_t>(nf);
  for (std::size_t i = 0; i < nf; ++i) {
    std::ifstream is(paths[i], std::ios::binary);
    if (!is) die("cannot open " + paths[i]);
    readHeaderFrom(is, set.headers[i]);
    if (set.headers[i].step != step)
      die(paths[i] + " holds step " + std::to_string(set.headers[i].step) +
          ", expected " + std::to_string(step));
    set.streams.push_back(std::move(is));
  }
  return set;
}

// compare two block records for the same block (geometric identity)
inline bool sameBlock(const BlockRecord& a, const BlockRecord& b) {
  return a.blockid == b.blockid && a.level == b.level && a.overlap == b.overlap &&
         a.dim == b.dim && a.mesh == b.mesh && a.n == b.n;
}

// merge the block tables of all rank files and align them with the current
// geometry by blockid:
//   - every local block must appear exactly once, in exactly one file, with an
//     identical geometric record
//   - the total block count must equal the local one (same domain, different
//     decomposition -- not a different domain)
// outputs fileToLive[fileIdx][blockIdx] = local block index (-1 if not ours)
std::string buildBlockMapping(const RankFileSet& files, const CheckpointHeader& live,
                              std::vector<std::vector<int>>& fileToLive) {
  // collect the blocks of every file
  std::vector<std::pair<std::size_t, std::size_t>> where;  // local block -> (fileIdx, blockIdx)
  std::vector<char> found(live.blocks.size(), 0);
  std::size_t totalFileBlocks = 0;
  fileToLive.assign(files.headers.size(), std::vector<int>{});
  for (std::size_t f = 0; f < files.headers.size(); ++f) {
    const auto& hdr = files.headers[f];
    fileToLive[f].assign(hdr.blocks.size(), -1);
    for (std::size_t b = 0; b < hdr.blocks.size(); ++b) {
      ++totalFileBlocks;
      for (std::size_t l = 0; l < live.blocks.size(); ++l) {
        if (sameBlock(hdr.blocks[b], live.blocks[l])) {
          if (found[l])
            return "block " + std::to_string(live.blocks[l].blockid) +
                   " appears in more than one rank's checkpoint file";
          found[l] = 1;
          fileToLive[f][b] = static_cast<int>(l);
          where.emplace_back(f, b);
          break;
        }
      }
    }
  }
  if (totalFileBlocks != live.blocks.size())
    return "the checkpoint holds " + std::to_string(totalFileBlocks) +
           " blocks in total, this run has " + std::to_string(live.blocks.size()) +
           " (a block cannot be split, so only a coarser-to-same decomposition "
           "can be resumed)";
  for (std::size_t l = 0; l < live.blocks.size(); ++l) {
    if (!found[l])
      return "block " + std::to_string(live.blocks[l].blockid) +
             " (mesh " + meshToString(live.blocks[l].mesh) + ") is missing from the " +
             "checkpoint; the decomposition or the geometry changed";
  }
  // the field metadata must match too (per field, matched by index)
  for (const auto& ff : files.headers[0].fields) {
    const FieldRecord* lf = nullptr;
    for (const auto& l : live.fields)
      if (l.index == ff.index) lf = &l;
    if (lf == nullptr) continue;  // extra field in the file: skip the payload, not an incompatibility
    if (lf->tag != ff.tag)
      return "field '" + ff.name + "' element type mismatch: file " +
             std::string(tagName(ff.tag)) + ", live " + std::string(tagName(lf->tag)) +
             " (rebuilt with different precision flags?)";
    if (lf->elem_dim != ff.elem_dim)
      return "field '" + ff.name + "' element dimension mismatch: file " +
             std::to_string(ff.elem_dim) + ", live " + std::to_string(lf->elem_dim);
    if (lf->name != ff.name)
      return "field " + std::to_string(ff.index) + " type mismatch: file '" + ff.name +
             "', live '" + lf->name + "'";
  }
  (void)where;
  return "";
}

// compare the file header with the current geometry item by item and return
// the first mismatch (empty string = compatible)
std::string explainMismatch(const CheckpointHeader& file, const CheckpointHeader& live) {
  if (file.dim != live.dim)
    return "dimension mismatch: file " + std::to_string(file.dim) + "D, live " +
           std::to_string(live.dim) + "D";
  if (file.blocks.size() != live.blocks.size())
    return "block count mismatch: file " + std::to_string(file.blocks.size()) +
           ", live " + std::to_string(live.blocks.size()) +
           " (restart requires the same decomposition and MPI size)";
  for (std::size_t i = 0; i < file.blocks.size(); ++i) {
    const auto& a = file.blocks[i];
    const auto& b = live.blocks[i];
    if (a.blockid != b.blockid)
      return "block " + std::to_string(i) + " id mismatch: file " +
             std::to_string(a.blockid) + ", live " + std::to_string(b.blockid);
    if (a.level != b.level)
      return "block " + std::to_string(a.blockid) + " level mismatch: file " +
             std::to_string(static_cast<int>(a.level)) + ", live " +
             std::to_string(static_cast<int>(b.level));
    if (a.overlap != b.overlap)
      return "block " + std::to_string(a.blockid) + " overlap mismatch: file " +
             std::to_string(a.overlap) + ", live " + std::to_string(b.overlap);
    if (a.n != b.n)
      return "block " + std::to_string(a.blockid) + " cell count mismatch: file " +
             std::to_string(a.n) + ", live " + std::to_string(b.n);
    if (a.mesh != b.mesh)
      return "block " + std::to_string(a.blockid) + " mesh mismatch: file " +
             meshToString(a.mesh) + ", live " + meshToString(b.mesh);
  }
  for (const auto& ff : file.fields) {
    const FieldRecord* lf = nullptr;
    for (const auto& l : live.fields)
      if (l.index == ff.index) lf = &l;
    if (lf == nullptr) continue;  // extra field in the file: skip the payload, not an incompatibility
    if (lf->tag != ff.tag)
      return "field '" + ff.name + "' element type mismatch: file " +
             std::string(tagName(ff.tag)) + ", live " + std::string(tagName(lf->tag)) +
             " (rebuilt with different precision flags?)";
    if (lf->elem_dim != ff.elem_dim)
      return "field '" + ff.name + "' element dimension mismatch: file " +
             std::to_string(ff.elem_dim) + ", live " + std::to_string(lf->elem_dim);
    if (lf->name != ff.name)
      return "field " + std::to_string(ff.index) + " type mismatch: file '" + ff.name +
             "', live '" + lf->name + "'";
    if (lf->counts != ff.counts)
      return "field '" + ff.name + "' array size mismatch (geometry changed?)";
  }
  return "";
}

}  // namespace detail

// ---------------------------------------------------------------------------
// public API
// ---------------------------------------------------------------------------
bool readHeader(const std::string& dir, std::int64_t step, CheckpointHeader& hdr) {
  if (step < 0) step = latestStep(dir);
  if (step < 0) {
    IF_MPI_RANK(0) { std::cerr << "[flbckpt] no checkpoint found in " << dir << std::endl; }
    return false;
  }
  const std::string path =
    dir + "/step_" + std::to_string(step) + "/rank" + std::to_string(mpi().getRank()) +
    ".ck";
  std::ifstream is(path, std::ios::binary);
  if (!is) {
    IF_MPI_RANK(0) { std::cerr << "[flbckpt] cannot open " << path << std::endl; }
    return false;
  }
  detail::readHeaderFrom(is, hdr);
  return true;
}

std::int64_t latestStep(const std::string& dir) {
  std::ifstream is(dir + "/latest");
  if (!is) return -1;
  std::int64_t step = -1;
  is >> step;
  return step;
}

template <typename... FMs>
void write(const std::string& dir, std::int64_t step, FMs&... fms) {
  const int rank = mpi().getRank();
  const std::string stepdir = dir + "/step_" + std::to_string(step);
  DirCreator::Create_Dir(dir);
  DirCreator::Create_Dir(stepdir);
  mpi().barrier();

  detail::MetaCollector mc;
  mc.mode = detail::MetaCollector::Mode::Write;
  int dummy[] = {0, (mc.run(fms), 0)...};
  (void)dummy;

  mc.live.step = step;
  mc.live.rank = rank;
  mc.live.nranks = mpi().getSize();

  // write a temporary file and rename it: a crash never leaves a partial
  // snapshot that looks complete
  const std::string tmp = stepdir + "/rank" + std::to_string(rank) + ".ck.tmp";
  const std::string dst = stepdir + "/rank" + std::to_string(rank) + ".ck";
  {
    std::ofstream os(tmp, std::ios::binary);
    if (!os) detail::die("cannot open " + tmp + " for writing");
    detail::writeHeader(os, mc.live);
    for (auto& w : mc.writers) w.second(os);
    os.flush();
    if (!os) detail::die("write failed for " + tmp);
  }
  if (std::rename(tmp.c_str(), dst.c_str()) != 0)
    detail::die("rename failed for " + tmp);
  mpi().barrier();

  // rank 0 writes "latest" last; its presence means this step's snapshot is
  // complete
  if (rank == 0) {
    const std::string latestTmp = dir + "/latest.tmp";
    {
      std::ofstream os(latestTmp, std::ios::trunc);
      os << step << "\n";
    }
    std::rename(latestTmp.c_str(), (dir + "/latest").c_str());
  }
  mpi().barrier();
}

template <typename... FMs>
std::int64_t read(const std::string& dir, std::int64_t step, FMs&... fms) {
  if (step < 0) step = latestStep(dir);
  if (step < 0) detail::die("no checkpoint found in " + dir);

  // ---- open **every** rank file of this step ----
  // difference from P1: P1 read only this rank's file (same rank count); here we
  // pick our blocks out of every rank file by blockid, so a different rank
  // count or a different decomposition can resume.
  detail::RankFileSet files = detail::openRankFiles(dir, step);

  detail::MetaCollector mc;
  mc.mode = detail::MetaCollector::Mode::Read;
  int dummy[] = {0, (mc.run(fms), 0)...};
  (void)dummy;

  // ---- validate and build the "file block -> local block" mapping ----
  std::vector<std::vector<int>> fileToLive;  // [fileIdx][blockIdx] = local block, -1 = not ours
  const std::string why = detail::buildBlockMapping(files, mc.live, fileToLive);
  if (!why.empty())
    detail::die("checkpoint " + dir + "/step_" + std::to_string(step) +
                " does not match this run: " + why);

  // ---- read the payload file by file, field by field, block by block ----
  for (std::size_t f = 0; f < files.streams.size(); ++f) {
    for (const auto& ff : files.headers[f].fields) {
      detail::MetaCollector::BlockReader reader = nullptr;
      for (auto& r : mc.readers)
        if (r.first == ff.index) reader = r.second;
      for (std::size_t b = 0; b < ff.counts.size(); ++b) {
        const int live = fileToLive[f][b];
        if (reader && live >= 0) {
          reader(files.streams[f], static_cast<std::size_t>(live), ff);
        } else {
          if (!reader) {
            IF_MPI_RANK(0) {
              std::cerr << "[flbckpt] warning: checkpoint field '" << ff.name
                        << "' is not present in this build; keeping its initial value"
                        << std::endl;
            }
          }
          detail::skipOneBlock(files.streams[f], ff, b);
        }
      }
    }
    if (!files.streams[f])
      detail::die("payload read failed for " + files.paths[f] + " (truncated file?)");
  }

#ifdef __CUDACC__
  for (auto& fn : mc.postRead) fn();
#endif
  mpi().barrier();
  return files.step;
}

template <typename... FMs>
bool checkCompatible(const std::string& dir, std::int64_t step, FMs&... fms) {
  if (step < 0) step = latestStep(dir);
  if (step < 0) return false;

  detail::MetaCollector mc;
  mc.mode = detail::MetaCollector::Mode::Check;
  int dummy[] = {0, (mc.run(fms), 0)...};
  (void)dummy;

  // open every rank file and run the same cross-decomposition validation
  // (without reading the payload)
  detail::RankFileSet files = detail::openRankFiles(dir, step, false);
  std::vector<std::vector<int>> fileToLive;
  const std::string why = detail::buildBlockMapping(files, mc.live, fileToLive);
  if (!why.empty()) {
    IF_MPI_RANK(0) {
      std::cerr << "[flbckpt] incompatible checkpoint " << dir << "/step_"
                << std::to_string(step) << ": " << why << std::endl;
    }
    return false;
  }
  return true;
}

}  // namespace flbckpt
