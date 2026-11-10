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

// aver_stencil.h
// One definition of "which fine cells feed a fine->coarse averaging", shared
// by the local (SharedComm / AddtoSharedAverComm) and the MPI
// (InitMPIAverComm) paths, in 2D and 3D.  Before this header the two paths --
// and the two MPI sides -- enumerated the averaging support independently, in
// two different coordinate frames, and for block pairs whose intersection is
// not a whole multiple of the level ratio the lengths did not agree at all.
//
// The rule:
//   the support is intersection(coarse FULL block, fine base block), counted
//   in COARSE cells.  The exchange writes the coarse block's overlap ring
//   facing the fine block, so the coarse side of the support includes that
//   ring; the fine side stops at the fine base (the cells that actually
//   carry data, not its halo).  For every coarse cell of the support the
//   stencil averages the 4 (2D) / 8 (3D) fine cells of its 2x2 (2x2x2) fine
//   neighbourhood; a fine coordinate that leaves the fine full block is
//   clamped into it ("the neighbourhood touches the boundary", the same
//   treatment the Intp exchange gets, see geometry/intp_stencil.h).
//
// Frames (chosen to reproduce the historical enumeration bit for bit on
// well-conditioned pairs):
//   coarse targets: the coarse block's FULL frame
//   fine sources:  the old path produced them from
//        coarsenedBlock().getRefinedCellIdx(), i.e. loc*2 inside a frame
//        anchored at FBase.min - 0.5*Fvox with mesh (base ext / (2*Fvox) + 1);
//        BuildAverStencil reproduces that frame directly (nC* below) without
//        materialising the coarsened block.
//
// Layouts written here (appended; caller reserves only):
//   local (SendRecvCells): per coarse cell [SrcNum fine ids..., 1 coarse id]
//        -- consumed by averCommunicate / BlockFieldManager::AverCommunicate
//        (stride AverSrcNum + 1, sources = neighbour/fine side).
//   mpi send (Sends): per coarse cell SrcNum fine ids, fine ids only
//        -- consumed by packAverSend (buffer = Cells.size() / SrcNum * q).
//   mpi recv (Recvs): per coarse cell 1 coarse id, coarse ids only
//        -- consumed by unpackCrossRecv (buffer = Cells.size() * q).

#pragma once

#include <cmath>
#include <array>
#include <cstddef>
#include <vector>

#include "geometry/basic_geometry.h"

namespace aver {

template <unsigned int D>
constexpr std::size_t SrcNum() {
  static_assert(D == 2 || D == 3, "Aver stencil only defined for 2D/3D");
  return D == 2 ? 4 : 8;
}

// emit the SrcN fine source ids of ONE coarse cell (c: coarsened-fine frame
// loc, refined to loc*2 + per-dimension offset; cc: the same cell in the
// coarse FULL frame), plus the coarse target for the layouts that want it
template <unsigned int D, typename CoarseIdx, typename RefinedStencil>
inline void EmitNeighbourhood(const Vector<int, D>& c, const Vector<int, D>& cc,
                              CoarseIdx coarseIdx, RefinedStencil refinedStencil,
                              std::vector<std::size_t>* SendRecvCells,
                              std::vector<std::size_t>* Sends, std::vector<std::size_t>* Recvs) {
  constexpr std::size_t SrcN = SrcNum<D>();
  const std::size_t cid = coarseIdx(cc);
  if (Recvs) Recvs->push_back(cid);
  const auto ids = refinedStencil(c);
  for (unsigned int s = 0; s < SrcN; ++s) {
    if (SendRecvCells) SendRecvCells->push_back(ids[s]);
    if (Sends) Sends->push_back(ids[s]);
  }
  // local layout: the coarse target follows its fine sources
  if (SendRecvCells) SendRecvCells->push_back(cid);
}

// CBase / FBase: the two base blocks of the coarse and the fine block of the
// pair (as stored in the geometry helper and in Block2D::getBaseBlock()).
// coarseOverlap / fineOverlap: the two blocks' overlaps; the full blocks are
// reconstructed the way Block2D does it.
// Exactly one of the output vectors may be null; the contents of the
// non-null ones are consistent across paths by construction:
//   Sends.size() == Recvs.size() * SrcNum
// and the local layout interleaves exactly these ids.
template <unsigned int D, typename T>
void BuildAverStencil(const BasicBlock<T, D>& CBase, const BasicBlock<T, D>& FBase,
                      int coarseOverlap, int fineOverlap,
                      std::vector<std::size_t>* SendRecvCells,
                      std::vector<std::size_t>* Sends, std::vector<std::size_t>* Recvs) {
  // the full (overlapped) blocks, same construction as Block2D(base, olap)
  const BasicBlock<T, D> CFull(CBase.getExtBlock(coarseOverlap));
  const BasicBlock<T, D> FFull(FBase.getExtBlock(fineOverlap));

  const T Fvox = FBase.getVoxelSize();
  const T Cvox = Fvox * T(2);

  // the support, counted in coarse cells over intersection(coarse FULL
  // block, fine base)
  const AABB<T, D> intsec = getIntersection(CFull, FBase);

  // coarsened-fine source frame: origin half a fine cell below the fine base
  // min, cell size 2*Fvox, mesh covering the base plus the anchor cell
  // (matches what getCoasenedBlock() produced: mesh = ext/Cvox + 1 with min
  // shifted by -0.5*Fvox)
  Vector<int, D> nCmesh{}, nCstart{}, nCnum{};
  for (unsigned int i = 0; i < D; ++i) {
    nCmesh[i] = static_cast<int>(std::round(FBase.getExtension()[i] / Cvox)) + 1;
    nCnum[i] = static_cast<int>(std::round(intsec.getExtension()[i] / Cvox));
    nCstart[i] =
      static_cast<int>(std::round((intsec.getMin()[i] - FBase.getMin()[i]) / Cvox + T(0.5)));
  }

  // coarse target index of the support start, in the coarse FULL frame
  Vector<int, D> startC{};
  for (unsigned int i = 0; i < D; ++i)
    startC[i] = static_cast<int>(std::round((intsec.getMin()[i] - CFull.getMin()[i]) / Cvox));

  auto coarseIdx = [&](const Vector<int, D>& loc) { return CFull.getIndex(loc); };
  // refined index in the coarsened-fine frame: loc*2 (matches
  // BasicBlock::getRefinedCellIdx on the coarsened block)
  auto refinedIdx = [&](const Vector<int, D>& loc) {
    if constexpr (D == 2) {
      const std::size_t m2 = static_cast<std::size_t>(nCmesh[0] * 2);
      return static_cast<std::size_t>(loc[0]) +
             static_cast<std::size_t>(loc[1]) * m2;
    } else {
      const std::size_t mx2 = static_cast<std::size_t>(nCmesh[0] * 2);
      const std::size_t XY = mx2 * static_cast<std::size_t>(nCmesh[1] * 2);
      return static_cast<std::size_t>(loc[0]) +
             static_cast<std::size_t>(loc[1]) * mx2 +
             static_cast<std::size_t>(loc[2]) * XY;
    }
  };
  // the neighbourhood offsets act on the REFINED mesh (loc*2 + off), matching
  // getRefinedCellIdx's {id0, id0+1, id0+row, id0+row+1} layout; coordinates
  // are clamped into the refined mesh (nCmesh*2) -- the neighbourhood touches
  // the boundary of the fine block's frame, the same treatment the Intp
  // stencil sources get
  auto refinedStencil = [&](const Vector<int, D>& c) {
    std::array<std::size_t, SrcNum<3>()> ids{};
    constexpr std::size_t SrcN = SrcNum<D>();
    const Vector<int, D> base = c * 2;
    for (unsigned int s = 0; s < SrcN; ++s) {
      Vector<int, D> p = base;
      for (unsigned int i = 0; i < D; ++i) {
        p[i] += (s >> i) & 1;
        const int m2 = nCmesh[i] * 2;
        p[i] = p[i] < 0 ? 0 : (p[i] >= m2 ? m2 - 1 : p[i]);
      }
      ids[s] = refinedIdx(p);
    }
    return ids;
  };
  (void)FFull;

  Vector<int, D> c{}, cc{};
  for (int iy = 0; iy < nCnum[1]; ++iy) {
    c[1] = iy + nCstart[1];
    cc[1] = iy + startC[1];
    for (int ix = 0; ix < nCnum[0]; ++ix) {
      c[0] = ix + nCstart[0];
      cc[0] = ix + startC[0];
      if constexpr (D == 3) {
        for (int iz = 0; iz < nCnum[2]; ++iz) {
          c[2] = iz + nCstart[2];
          cc[2] = iz + startC[2];
          EmitNeighbourhood<D>(c, cc, coarseIdx, refinedStencil, SendRecvCells, Sends, Recvs);
        }
      } else {
        EmitNeighbourhood<D>(c, cc, coarseIdx, refinedStencil, SendRecvCells, Sends, Recvs);
      }
    }
  }
  (void)startC;
}

}  // namespace aver
