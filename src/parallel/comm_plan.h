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

// comm_plan.h

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include "parallel/link_key.h"

namespace flbcomm {

// Which communicator set an exchange uses.  Partial is the direction-filtered
// DirSends/DirRecvs pair, Full is MPIComm and All is the wider AllMPIComm.
// The three sets have their own cell lists and message sizes; the hidden step
// only needs to know which plan and which pack/unpack to use.
enum class CommMode { Partial, Full, All };

// Cross-level exchanges.  Aver (fine -> coarse) and Intp (coarse -> fine) have
// their own cell lists and message layout; the All variants use the wider
// overlap set.  Each is a separate call, so it only needs one tag slot.
enum class CrossMode { Aver, Intp, AllAver, AllIntp };

inline constexpr bool crossIsAll(CrossMode mode) {
  return mode == CrossMode::AllAver || mode == CrossMode::AllIntp;
}
inline constexpr bool crossIsAver(CrossMode mode) {
  return mode == CrossMode::Aver || mode == CrossMode::AllAver;
}

// One message per peer rank is posted per plan, so the MPI tag only has to
// distinguish concurrently pending exchanges, not individual links.  The links
// inside a plan are ordered by CommRef::Order.  Keeping the tags small also
// removes the old need to derive them from a global block/link counter, which
// silently approached MPI_TAG_UB on large decompositions.
inline constexpr int DirTag = 1;
inline constexpr int FullTag = 2;
inline constexpr int AllTag = 3;
inline constexpr int AverTag = 4;
inline constexpr int IntpTag = 5;
inline constexpr int AllAverTag = 6;
inline constexpr int AllIntpTag = 7;

inline constexpr int crossTag(CrossMode mode) {
  return mode == CrossMode::Aver      ? AverTag
         : mode == CrossMode::Intp    ? IntpTag
         : mode == CrossMode::AllAver ? AllAverTag
                                      : AllIntpTag;
}

// Field plans start here, above every lattice tag.  A field exchange can be in
// flight while a lattice exchange is pending (the hidden step), and MPI matches
// on (communicator, source, tag), so the ranges must not overlap.
inline constexpr int FieldTagBase = 1 << 14;
// Each field manager reserves six slots (normal, all, aver, intp, all-aver,
// all-intp).  The bases are handed out in construction order, which every rank
// must follow identically: the application builds the same field managers
// everywhere.
inline constexpr int FieldTagStride = 6;
inline constexpr int FieldTagSpace = 1 << 12;

inline constexpr int fieldTagSlot(CrossMode mode) {
  return mode == CrossMode::Aver      ? 2
         : mode == CrossMode::Intp    ? 3
         : mode == CrossMode::AllAver ? 4
                                      : 5;
}

inline int allocateFieldTagBase() {
  static int next = FieldTagBase;
  const int base = next;
  next += FieldTagStride;
  if (next > FieldTagBase + FieldTagSpace) {
    std::cerr << "[flbcomm] field tag range exhausted (FieldTagSpace = "
              << FieldTagSpace << ").\n";
    std::exit(1);
  }
  return base;
}

// Shape of an Aver/Intp communication list: each entry is `SrcNum` source
// cells followed by one destination cell, so a pair is `Stride` long.  The
// lists are built by AddtoSharedAverComm() / AddtoSharedIntpComm().
inline constexpr std::size_t AverSrcNum(std::size_t dim) { return dim == 2 ? 4 : 8; }
inline constexpr std::size_t AverPairStride(std::size_t dim) { return dim == 2 ? 5 : 9; }

// Elements in one cross-level message.  A send packs q (or a field's array_dim)
// values per AverSrcNum source cells; a receive takes every listed cell.
inline std::size_t crossSendElements(std::size_t cells, std::size_t dim,
                                     std::size_t q) {
  return cells / AverSrcNum(dim) * q;
}
inline std::size_t crossRecvElements(std::size_t cells, std::size_t q) {
  return cells * q;
}

// Level filter for the sub-cycled exchange.  A block whose `deLevel`
// (maxLevel - level) satisfies 2^deLevel | count takes part, which is the same
// as deLevel <= ctz(count).  The active set therefore depends only on the
// trailing zeros of count, capped at maxLevel -- at most maxLevel + 1 distinct
// sets, so a caller can cache one plan per depth instead of rebuilding.
inline std::size_t activeDepth(std::int64_t count, std::size_t maxLevel) {
  if (count <= 0) return maxLevel;  // 0 % k == 0 for every k
  std::size_t v = 0;
  while (v < maxLevel && ((count >> v) & 1) == 0) ++v;
  return v;
}

// One logical link (a single DistributedComm of a single block) inside a
// coalesced message.  Block and Comm index into the owner's block vector and
// that block's comm list; Offset is where the link starts inside the group
// buffer.
struct CommRef {
  std::size_t Block{};
  std::size_t Comm{};
  std::size_t Offset{};
  int Rank{};  // target rank for sends, source rank for receives
  int Tag{};
  // Total order inside the group.  The builder picks a key that both sides can
  // derive from the geometry alone; the refs of a group are sorted by it.
  LinkKey Order{};
  // Elements this ref contributes to its group buffer.  The engine uses it to
  // compute offsets; managers no longer pass a size lambda to buildGroups().
  std::size_t Size{};
};

// One MPI message: every link to or from the same peer, packed back to back.
// The sender and the receiver build their plans independently, so the ref
// order must be derivable from the geometry alone (see the builders in
// block_lattice.hh).
template <typename buffer_type>
struct CommGroup {
  int TargetRank{};
  int Tag{};  // smallest tag in the group; both sides pick the same one
  std::vector<CommRef> Refs;
  buffer_type Buffer;  // owned here so it outlives the non-blocking call
};

// Send and receive plans for one communication path.  The buffers live in the
// plan, so a step only resizes them instead of allocating a buffer per block
// per communicator.
template <typename buffer_type>
struct CommPlan {
  std::vector<CommGroup<buffer_type>> Send;
  std::vector<CommGroup<buffer_type>> Recv;
  bool Ready = false;
  // Geometry generation this plan was built against.  The stamps are per plan,
  // not per manager: one manager owns plans for several depths and exchange
  // kinds, and a single rebuild must not vouch for the plans it did not touch.
  std::uint64_t Generation{};
};

// Sort the refs by (peer rank, packing order) and cut them into one group per
// peer, packing the refs of a group back to back.  Each ref carries its own
// Size; the group tag is the smallest tag in it, which is what both sides agree
// on without negotiating.
template <typename buffer_type>
void buildGroups(std::vector<CommRef>& refs,
                 std::vector<CommGroup<buffer_type>>& groups) {
  std::sort(refs.begin(), refs.end(), [](const CommRef& a, const CommRef& b) {
    return a.Rank != b.Rank ? a.Rank < b.Rank : a.Order < b.Order;
  });
  groups.clear();
  for (std::size_t i = 0; i < refs.size();) {
    std::size_t j = i + 1;
    while (j < refs.size() && refs[j].Rank == refs[i].Rank) ++j;
    CommGroup<buffer_type> group;
    group.TargetRank = refs[i].Rank;
    std::size_t offset = 0;
    group.Refs.reserve(j - i);
    for (std::size_t k = i; k < j; ++k) {
      group.Tag = k == i ? refs[k].Tag : std::min(group.Tag, refs[k].Tag);
      refs[k].Offset = offset;
      offset += refs[k].Size;
      group.Refs.push_back(refs[k]);
    }
    group.Buffer.resize(offset);
    groups.push_back(std::move(group));
    i = j;
  }
}

// Per-side variant for exchanges whose send and receive sets are filtered
// differently (the cross-level Aver/Intp sub-cycling).  `active(block, send)`
// decides whether that side of the block contributes refs.  `generation` is
// the geometry generation the caller is building against and is stored on the
// plan, so every plan carries the stamp of its own build.
template <typename buffer_type, typename ActiveFn, typename LinkFn>
void buildPlan(std::size_t blockCount, std::uint64_t generation, ActiveFn&& active,
               LinkFn&& linksForBlock, CommPlan<buffer_type>& plan) {
  std::vector<CommRef> sendRefs;
  std::vector<CommRef> recvRefs;
  for (std::size_t bi = 0; bi < blockCount; ++bi) {
    if (active(bi, true)) linksForBlock(bi, true, sendRefs);
    if (active(bi, false)) linksForBlock(bi, false, recvRefs);
  }
  buildGroups(sendRefs, plan.Send);
  buildGroups(recvRefs, plan.Recv);
  plan.Generation = generation;
  plan.Ready = true;
}

// Same-level variant: a block takes part on both sides when it is within
// `depth` of the finest level.  depth == maxLevel is the unfiltered exchange.
template <typename buffer_type, typename LevelFn, typename LinkFn>
void buildPlan(std::size_t blockCount, std::size_t maxLevel, std::size_t depth,
               std::uint64_t generation, LevelFn&& levelOf, LinkFn&& linksForBlock,
               CommPlan<buffer_type>& plan) {
  buildPlan<buffer_type>(
    blockCount, generation,
    [&](std::size_t bi, bool) {
      const std::size_t level = static_cast<std::size_t>(levelOf(bi));
      return level <= maxLevel && maxLevel - level <= depth;
    },
    linksForBlock, plan);
}

// A plan indexes cell lists built by the geometry, so a rebuild of that
// geometry leaves it pointing at cells that no longer hold what it was built
// against.  Every consumer compares the plan's own stamp with the current
// generation before using it; the manager-wide "last built" value this
// replaces let any later build revalidate the older plans.
inline void checkPlanStamp(std::uint64_t planGeneration,
                           std::uint64_t currentGeneration, const char* owner) {
  if (planGeneration == currentGeneration) return;
  std::cerr << "[flbcomm] " << owner << ": communication plan is stale; built for "
               "geometry generation "
            << planGeneration << ", the geometry is at " << currentGeneration
            << ".  Rebuild the communicators after the geometry changed.\n";
#ifdef MPI_ENABLED
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  std::exit(1);
}

// Guard for the tag namespace.  Every exchange uses a fixed slot or a
// per-manager range now, so only the reserved ranges have to fit under
// MPI_TAG_UB; block ids do not take part.
inline void checkTagBudget(int tagUpperBound) {
  const int highestReserved = FieldTagBase + FieldTagSpace - 1;
  if (tagUpperBound < highestReserved) {
    std::cerr << "[flbcomm] MPI_TAG_UB is " << tagUpperBound
              << ", below the reserved tag range [" << FieldTagBase << ", "
              << highestReserved << "].\n";
    std::exit(1);
  }
}

// Guard for a cell index that a communication list claims to address.  The
// Aver/Intp lists are built by the geometry (`AddtoSharedAverComm` /
// `AddtoSharedIntpComm`) and the consumer writes through them unconditionally,
// so a list that reaches past the end of its block used to corrupt the heap of
// whatever allocation happened to follow.  Report it instead: the index is
// wrong, not the data, and silently skipping it would hide a real geometry bug.
inline void checkCellIndex(std::size_t id, std::size_t n, int blockId, int peerBlockId,
                           const char* what) {
  if (id < n) return;
  std::cerr << "[flbcomm] " << what << ": cell index " << id << " is outside block "
            << blockId << " (N = " << n << "), neighbour block " << peerBlockId
            << ".  The Aver/Intp communication list of the geometry is inconsistent.\n";
  std::exit(1);
}

}  // namespace flbcomm
