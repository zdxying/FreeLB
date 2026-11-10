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
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 *
 */

// intp_stencil.h
// One definition of "which coarse cells take part in a coarse->fine
// interpolation", shared by the local (SharedComm / AddtoSharedIntpComm) and
// the MPI (InitMPIIntpComm) paths, in 2D and 3D.  Before this header the two
// paths enumerated the stencils independently and the two MPI sides disagreed
// on the set -- send length != recv length, the exchange deadlocked.
//
// The rule, matching the rows of getIntpWeight<D>():
//   the support is intersection(fine full block, coarse base block), counted
//   in coarse cells; every coarse unit expands to 4 (2D) / 8 (3D) stencils,
//   each reading a 2x2 (2x2x2) coarse neighbourhood and writing one fine
//   cell.  A coordinate that leaves its block is clamped into the block's
//   FULL extent (halo included) -- "the stencil touches the boundary", the
//   same treatment the Aver exchange gets from the geometry.  There is
//   deliberately no source outside a block: the MPI sides must enumerate the
//   identical set, and they both run this function on the same two base
//   blocks, so send == recv * SrcNum holds by construction.
//
// Layouts written here (appended; caller reserves only):
//   local (SendRecvCells): per stencil [SrcNum coarse ids..., 1 fine id],
//        NumStencils per coarse unit  -- consumed by popIntpCommunicate.
//   mpi send (Sends): per coarse unit NumStencils * SrcNum coarse ids,
//        stencil-major  -- consumed by packIntpSend / popIntp.
//   mpi recv (Recvs): per coarse unit NumStencils fine ids, one per stencil
//        in the same shift order  -- consumed by unpackCrossRecv.

#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "geometry/basic_geometry.h"

namespace intp {

template <unsigned int D>
constexpr std::size_t SrcNum() {
  static_assert(D == 2 || D == 3, "Intp stencil only defined for 2D/3D");
  return D == 2 ? 4 : 8;
}

// emit the NumStencils stencils of ONE coarse unit: for shift s the source
// neighbourhood sits at c0 + shift(s) .. + 1, the target fine cell at
// f + shift(s); shift order matches the rows of getIntpWeight<D>()
template <unsigned int D, typename T, typename FineIdx, typename CoarseIdx>
inline void EmitStencil(const Vector<int, D>& f, const Vector<int, D>& c0,
                        const BasicBlock<T, D>& FFull, const BasicBlock<T, D>& CFull,
                        FineIdx fineIdx, CoarseIdx coarseIdx,
                        std::vector<std::size_t>* SendRecvCells,
                        std::vector<std::size_t>* Sends, std::vector<std::size_t>* Recvs) {
  constexpr std::size_t SrcN = SrcNum<D>();
  constexpr int NumStencils = 1 << D;

  for (int shift = 0; shift < NumStencils; ++shift) {
    Vector<int, D> fshift{}, cshift{};
    for (unsigned int i = 0; i < D; ++i) {
      fshift[i] = (shift >> i) & 1;
      cshift[i] = fshift[i];
    }
    // fine target of this stencil: one entry per stencil on the MPI side
    if (Recvs) Recvs->push_back(fineIdx(f + fshift));
    // SrcN coarse sources of this stencil
    for (unsigned int s = 0; s < SrcN; ++s) {
      Vector<int, D> soff{};
      for (unsigned int i = 0; i < D; ++i) soff[i] = (s >> i) & 1;
      const std::size_t cid = coarseIdx(c0 + cshift + soff);
      if (SendRecvCells) SendRecvCells->push_back(cid);
      if (Sends) Sends->push_back(cid);
    }
    // local layout: the fine target follows its stencil's sources
    if (SendRecvCells) SendRecvCells->push_back(fineIdx(f + fshift));
  }
}

// FBase / CBase: the two base blocks of the fine and the coarse block of the
// pair (as stored in the geometry helper and in Block2D::getBaseBlock()).
// overlap: the block overlap; the full blocks are reconstructed from the
// bases exactly the way Block2D does it.
// Exactly one of the output vectors may be null; the contents of the
// non-null ones are consistent across paths by construction.
template <unsigned int D, typename T>
void BuildStencil(const BasicBlock<T, D>& FBase, const BasicBlock<T, D>& CBase,
                  int fineOverlap, int coarseOverlap,
                  std::vector<std::size_t>* SendRecvCells,
                  std::vector<std::size_t>* Sends, std::vector<std::size_t>* Recvs) {
  // the full (overlapped) blocks, same construction as Block2D(base, olap);
  // the two blocks' overlaps can differ (the level rule gives 2 above level 0)
  const BasicBlock<T, D> FFull(FBase.getExtBlock(fineOverlap));
  const BasicBlock<T, D> CFull(CBase.getExtBlock(coarseOverlap));

  const T Fvox = FFull.getVoxelSize();
  const T Cvox = CFull.getVoxelSize();

  // the support, counted in coarse cells over intersection(fine full, coarse base)
  const AABB<T, D> intsec = getIntersection(FFull, CBase);
  Vector<int, D> CN{}, startF{}, startC{};
  for (unsigned int i = 0; i < D; ++i) {
    CN[i] = static_cast<int>(std::round(intsec.getExtension()[i] / Cvox));
    startF[i] = static_cast<int>(std::round((intsec.getMin()[i] - FFull.getMin()[i]) / Fvox));
    // shift 1 coarse voxel to the left bottom for the interpolation templates
    startC[i] = static_cast<int>(std::round((intsec.getMin()[i] - CFull.getMin()[i]) / Cvox)) - 1;
  }

  // clamp-in-block index; the clamp bound is the FULL block (halo included)
  auto fineIdx = [&](const Vector<int, D>& loc) {
    Vector<int, D> l;
    for (unsigned int i = 0; i < D; ++i)
      l[i] = loc[i] < 0 ? 0 : (loc[i] >= FFull.getMesh()[i] ? FFull.getMesh()[i] - 1 : loc[i]);
    return FFull.getIndex(l);
  };
  auto coarseIdx = [&](const Vector<int, D>& loc) {
    Vector<int, D> l;
    for (unsigned int i = 0; i < D; ++i)
      l[i] = loc[i] < 0 ? 0 : (loc[i] >= CFull.getMesh()[i] ? CFull.getMesh()[i] - 1 : loc[i]);
    return CFull.getIndex(l);
  };

  Vector<int, D> f{}, c0{};
  for (int iy = 0; iy < CN[1]; ++iy) {
    f[1] = iy * 2 + startF[1];
    c0[1] = iy + startC[1];
    for (int ix = 0; ix < CN[0]; ++ix) {
      f[0] = ix * 2 + startF[0];
      c0[0] = ix + startC[0];
      if constexpr (D == 3) {
        for (int iz = 0; iz < CN[2]; ++iz) {
          f[2] = iz * 2 + startF[2];
          c0[2] = iz + startC[2];
          EmitStencil<D>(f, c0, FFull, CFull, fineIdx, coarseIdx, SendRecvCells, Sends, Recvs);
        }
      } else {
        EmitStencil<D>(f, c0, FFull, CFull, fineIdx, coarseIdx, SendRecvCells, Sends, Recvs);
      }
    }
  }
}

}  // namespace intp
