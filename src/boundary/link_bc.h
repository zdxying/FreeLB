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
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// link_bc.h
// Static interpolated bounce-back on the narrow-band surface links.  The
// manager builds a block-local link table from VoxelGeometry and the flag
// field, then applies the Bouzidi two-point / one-sided formulas after Stream()
// and the cell dynamics.  Host only: the device path still uses the flag-based
// bounce-back task.

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

#include "boundary/surface_link.h"
#include "data_struct/Vector.h"
#include "geometry/voxel_geometry.h"
#include "lbm/lattice_set.h"
#include "utils/util.h"

namespace flbbc {

// Rigid rotation of a wall entity about an axis.  The link builder evaluates
// u_w = Omega * (Axis x (x_wall - Center)) once per link, because the wall is
// rigid: the velocity does not change with time.  A constant per-entity table
// (WallVelocity) cannot express a rotating surface, whose speed grows with r.
template <typename T>
struct WallMotion {
  Vector<T, 3> Center{};
  Vector<T, 3> Axis{T(0), T(0), T(1)};
  T Omega{};
};

template <typename BLOCKLATTICEMANAGER, typename BLOCKFIELDMANAGER>
class InterpolatedBBManager {
 public:
  using LatSet = typename BLOCKLATTICEMANAGER::LatticeSet;
  using T = typename BLOCKLATTICEMANAGER::FloatType;
  using CELL = typename BLOCKLATTICEMANAGER::CellType;
  using GenericRho = typename CELL::GenericRho;

  // periodicAxes: bit 0 = x, 1 = y, 2 = z.  A periodic cell keeps the links
  // that do not leave the domain along one of these axes; links across the
  // periodic face are left to the copy.  0 keeps the old whole-cell skip.
  InterpolatedBBManager(
    BLOCKLATTICEMANAGER& lat, BLOCKFIELDMANAGER& flagFM,
    const flbgeo::VoxelGeometry<T>& record,
    std::uint8_t voidflag = std::uint8_t{1},
    std::uint8_t periodicFlag = std::uint8_t{0},
    const std::unordered_map<std::uint16_t, Vector<T, 3>>* wallVelocities = nullptr,
    const std::unordered_map<std::uint16_t, WallMotion<T>>* wallMotions = nullptr,
    std::uint8_t periodicAxes = std::uint8_t{0})
      : _lat(lat), _flagFM(flagFM), _record(record), _voidflag(voidflag),
        _periodicFlag(periodicFlag), _wallVelocities(wallVelocities),
        _wallMotions(wallMotions), _periodicAxes(periodicAxes) {
    Build();
  }

  void Apply() {
    for (std::size_t ib = 0; ib < _links.size(); ++ib) {
      auto& blocklat = _lat.getBlockLat(static_cast<int>(ib));
      const auto& links = _links[ib];
      // Links of one cell are contiguous (Build emits them cell by cell).
      std::size_t i = 0;
      while (i < links.size()) {
        const std::size_t cellId = links[i].CellId;
        std::size_t end = i;
        while (end < links.size() && links[end].CellId == cellId) ++end;
        CELL cell(cellId, blocklat);
        T massError = T(0);
        for (std::size_t j = i; j < end; ++j) {
          const SurfaceLink<T>& link = links[j];
          const unsigned int k = link.Direction;
          const unsigned int opp = latset::opp<LatSet>(k);
          const T q = link.Fraction;
          T wall = T(0);
          if (link.WallVelocity[0] != T(0) || link.WallVelocity[1] != T(0) ||
              link.WallVelocity[2] != T(0)) {
            const auto copp = latset::c<LatSet>(opp);
            const T uc = link.WallVelocity[0] * T(copp[0]) +
                         link.WallVelocity[1] * T(copp[1]) +
                         link.WallVelocity[2] * T(copp[2]);
            // The q=1/2 form is 2 w rho (c_opp.u_w) / cs^2, the same term as
            // movingwall_bounceback.  The one-sided branch divides by q; the
            // two-point branch keeps the factor 2.
            const T scale = (q < T(0.5) && !link.Thin) ? T(2) : (T(1) / q);
            wall = scale * LatSet::InvCs2 * latset::w<LatSet>(opp) *
                   cell.template get<GenericRho>() * uc;
          }
          const T old = cell[opp];
          const T fk = cell.getPrevious(k);
          if (q < T(0.5) && !link.Thin) {
            CELL back = cell.getNeighbor(opp);
            const T fk2 = back.getPrevious(k);
            cell[opp] = T(2) * q * fk + (T(1) - T(2) * q) * fk2 + wall;
          } else {
            const T fopp = cell.getPrevious(opp);
            cell[opp] = fk / (T(2) * q) +
                        (T(2) * q - T(1)) / (T(2) * q) * fopp + wall;
          }
          massError += cell[opp] - old;
        }
        // The interpolated reflection replaces a population and, for a moving
        // wall, adds a term.  On a complete flat wall the per-cell change
        // cancels between the symmetric links; a strip truncated by a periodic
        // face leaves a net mass change that accumulates over the run.  w_k is
        // the zero-momentum mode (sum w = 1, sum w c = 0), so removing the
        // change from every direction conserves mass without touching
        // momentum.
        if (massError != T(0)) {
          for (unsigned int k = 0; k < LatSet::q; ++k) {
            cell[k] -= massError * latset::w<LatSet>(k);
          }
        }
        i = end;
      }
    }
  }

  std::size_t getLinkCount() const {
    std::size_t n = 0;
    for (const auto& links : _links) n += links.size();
    return n;
  }

  // Why a link fell back to the half-way position: no band sample at the
  // fluid cell, the band normal does not cross the link, or q left (0, 1].
  // Fallback links are valid but less accurate, so the split guides the
  // curved-wall work.
  std::size_t getFallbackNoSample() const { return _fallbackNoSample; }
  std::size_t getFallbackNoCross() const { return _fallbackNoCross; }
  std::size_t getFallbackOutOfRange() const { return _fallbackOutOfRange; }

  const std::vector<SurfaceLink<T>>& getLinks(int block) const {
    return _links[static_cast<std::size_t>(block)];
  }

 private:
  void Build() {
    const int blockNum = static_cast<int>(_lat.getBlockLats().size());
    _links.assign(static_cast<std::size_t>(blockNum), {});
    const T voxel = _record.getVoxelSize();
    const T eps = T(16) * std::numeric_limits<T>::epsilon();
    for (int ib = 0; ib < blockNum; ++ib) {
      auto& blocklat = _lat.getBlockLat(ib);
      const auto& block = blocklat.getBlock();
      const auto& flag = _flagFM.getBlockField(ib).getField(0);
      const int overlap = blocklat.getOverlap();
      const Vector<int, 3> mesh = block.getBaseBlock().getMesh();
      const int nx = blocklat.getNx();
      const int ny = blocklat.getNy();
      const int nz = blocklat.getNz();
      for (int z = overlap; z < overlap + mesh[2]; ++z) {
        for (int y = overlap; y < overlap + mesh[1]; ++y) {
          for (int x = overlap; x < overlap + mesh[0]; ++x) {
            const Vector<int, 3> loc{x, y, z};
            const std::size_t id = block.getIndex(loc);
            if (flag[id] == _voidflag) continue;
            // A periodic face is refreshed by the periodic manager, so a link
            // that would leave the domain across that face is skipped.  The
            // other directions stay: a wall that touches a periodic face still
            // needs its wall link on that face column, and the copy alone
            // cannot supply a link that was never built.
            const bool periodic =
              _periodicFlag != 0 && util::isFlag(flag[id], _periodicFlag);
            if (periodic && _periodicAxes == 0) continue;  // legacy whole-cell skip
            const Vector<T, 3> pos = block.getVoxel(loc);
            const std::size_t gid = _record.findCellId(pos);
            const flbgeo::SurfaceSample<T>* sample = nullptr;
            if (gid != flbgeo::kInvalidCell &&
                _record.hasSurfaceSample(gid)) {
              sample = &_record.getSurfaceSample(gid);
            }
            // A closed STL may enclose the fluid (geo_stl) or a solid (a
            // sphere): signedDistance is negative inside the mesh and
            // evalSurfaceNormal points away from the fluid in the first case.
            // Normalise both to a positive wall distance and a normal that
            // points from the wall into the fluid.
            T wallDistance = T(0);
            Vector<T, 3> nFluid{};
            if (sample != nullptr) {
              wallDistance = std::abs(sample->Distance);
              if (!(wallDistance > eps)) {
                sample = nullptr;
              } else {
                const bool fluidInside = sample->Distance < T(0);
                nFluid =
                  fluidInside
                    ? Vector<T, 3>{-sample->Normal[0], -sample->Normal[1],
                                   -sample->Normal[2]}
                    : sample->Normal;
              }
            }
            for (unsigned int k = 1; k < LatSet::q; ++k) {
              const auto c = latset::c<LatSet>(k);
              const Vector<int, 3> nloc{x + c[0], y + c[1], z + c[2]};
              if (nloc[0] < 0 || nloc[0] >= nx || nloc[1] < 0 ||
                  nloc[1] >= ny || nloc[2] < 0 || nloc[2] >= nz) {
                continue;
              }
              const std::size_t nid = block.getIndex(nloc);
              if (flag[nid] != _voidflag) continue;
              if (periodic) {
                const Vector<T, 3> npos = block.getVoxel(nloc);
                const auto& dmin = _record.getDomain().getMin();
                const auto& dmax = _record.getDomain().getMax();
                bool leavesPeriodic = false;
                for (int axis = 0; axis < 3; ++axis) {
                  if ((_periodicAxes & std::uint8_t(1u << axis)) == 0) continue;
                  if (npos[axis] < dmin[axis] || npos[axis] >= dmax[axis]) {
                    leavesPeriodic = true;
                    break;
                  }
                }
                // This direction crosses a periodic face; the copy owns it.
                if (leavesPeriodic) continue;
              }
              SurfaceLink<T> link;
              link.CellId = id;
              link.Direction = k;
              link.Entity = sample != nullptr ? sample->Entity : kNoEntity;
              T q = T(0.5);
              bool fallback = true;
              if (sample != nullptr) {
                const T dot = nFluid[0] * T(c[0]) + nFluid[1] * T(c[1]) +
                              nFluid[2] * T(c[2]);
                if (dot < -eps) {
                  q = wallDistance / (voxel * (-dot));
                  if (q > T(1)) {
                    // The void neighbour proves the wall is before the link
                    // end; clamping is closer than collapsing to the half-way
                    // position.  q <= 0 still means the sample is unusable.
                    q = T(1);
                    fallback = false;
                    ++_fallbackOutOfRange;
                  } else if (!(q > T(0))) {
                    fallback = true;
                    ++_fallbackOutOfRange;
                  } else {
                    fallback = false;
                  }
                } else {
                  ++_fallbackNoCross;
                }
              } else {
                ++_fallbackNoSample;
              }
              if (fallback) q = T(0.5);
              link.Fraction = q;
              link.Fallback = fallback;
              if (sample != nullptr) {
                bool setByMotion = false;
                if (_wallMotions != nullptr) {
                  auto it = _wallMotions->find(sample->Entity);
                  if (it != _wallMotions->end()) {
                    const Vector<T, 3> wallPoint = pos - wallDistance * nFluid;
                    link.WallVelocity =
                      it->second.Omega * CrossProduct3(it->second.Axis,
                                                       wallPoint - it->second.Center);
                    setByMotion = true;
                  }
                }
                if (!setByMotion && _wallVelocities != nullptr) {
                  auto it = _wallVelocities->find(sample->Entity);
                  if (it != _wallVelocities->end()) {
                    link.WallVelocity = it->second;
                  }
                }
              }
              if (!fallback && q < T(0.5)) {
                const Vector<int, 3> back{x - c[0], y - c[1], z - c[2]};
                const bool inside =
                  back[0] >= 0 && back[0] < nx && back[1] >= 0 &&
                  back[1] < ny && back[2] >= 0 && back[2] < nz;
                if (!inside || flag[block.getIndex(back)] == _voidflag) {
                  link.Thin = true;
                  link.Fraction = T(0.5);
                }
              }
              _links[static_cast<std::size_t>(ib)].push_back(link);
            }
          }
        }
      }
    }
  }

  BLOCKLATTICEMANAGER& _lat;
  BLOCKFIELDMANAGER& _flagFM;
  const flbgeo::VoxelGeometry<T>& _record;
  std::uint8_t _voidflag;
  std::uint8_t _periodicFlag;
  const std::unordered_map<std::uint16_t, Vector<T, 3>>* _wallVelocities;
  const std::unordered_map<std::uint16_t, WallMotion<T>>* _wallMotions;
  std::uint8_t _periodicAxes = 0;
  std::size_t _fallbackNoSample = 0;
  std::size_t _fallbackNoCross = 0;
  std::size_t _fallbackOutOfRange = 0;
  std::vector<std::vector<SurfaceLink<T>>> _links;
};

}  // namespace flbbc
