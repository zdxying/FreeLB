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

// link_bc_curved.cpp
// B4 curved-wall accuracy: Taylor-Couette flow between two concentric
// cylinders, periodic along z so the analytic profile has no end effects.  The
// inner cylinder rotates as a rigid body (entity 2), the outer shell is at rest
// (entity 3).  For a Stokes flow the azimuthal velocity is
//   v_theta(r) = A r + B / r,
//   A = -Omega ri^2 / (ro^2 - ri^2),  B = Omega ri^2 ro^2 / (ro^2 - ri^2),
// which is zero at ro and Omega * ri at ri.  The numerical profile at the
// mid-plane is compared with it cell by cell.
//
// Both meshes are written as ASCII STL next to the binary and removed after the
// geometry record is built; the test needs no committed fixture.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "freelb.h"
#include "freelb.hh"
#include "boundary/link_bc.h"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Thread_Num = 1;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);
constexpr std::uint8_t PeriodicFlag = std::uint8_t(4);
constexpr int Nx = 24;
constexpr int Ny = 24;
constexpr int Nz = 6;
constexpr T Center = T(12);
constexpr T InnerRadius = T(5);
constexpr T OuterRadius = T(10);
constexpr T ShellOuterRadius = T(18);
constexpr T Omega = T(0.0005);
constexpr int Steps = 6000;
constexpr int ThetaSegments = 72;
constexpr T PiConst = T(3.14159265358979323846);

using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL = Cell<T, LatSet, FIELDS>;
using BulkTask =
  tmp::Key_TypePair<BulkFlag,
                    collision::BGK<moment::rhoU<CELL>,
                                   equilibrium::SecondOrder<CELL>>>;
using PeriodicBulkTask =
  tmp::Key_TypePair<BulkFlag | PeriodicFlag,
                    collision::BGK<moment::rhoU<CELL>,
                                   equilibrium::SecondOrder<CELL>>>;
using TaskCollection = tmp::TupleWrapper<BulkTask, PeriodicBulkTask>;
using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;
using RhoUTask = tmp::Key_TypePair<BulkFlag, moment::rhoU<CELL, true>>;
using PeriodicRhoUTask =
  tmp::Key_TypePair<BulkFlag | PeriodicFlag, moment::rhoU<CELL, true>>;
using RhoUTaskCollection = tmp::TupleWrapper<RhoUTask, PeriodicRhoUTask>;
using TaskSelectorRhoU =
  tmp::TaskSelector<RhoUTaskCollection, std::uint8_t, CELL>;

using Tri = std::array<Vector<T, 3>, 3>;

Vector<T, 3> polar(T radius, T angle, T z) {
  return Vector<T, 3>{Center + radius * std::cos(angle),
                      Center + radius * std::sin(angle), z};
}

// Solid cylinder: side normals point away from the axis, caps up/down.
void addCylinder(std::vector<Tri>& tris, T radius, T z0, T z1) {
  for (int i = 0; i < ThetaSegments; ++i) {
    const T a = T(2) * PiConst * T(i) / T(ThetaSegments);
    const T b = T(2) * PiConst * T(i + 1) / T(ThetaSegments);
    const Vector<T, 3> a0 = polar(radius, a, z0);
    const Vector<T, 3> b0 = polar(radius, b, z0);
    const Vector<T, 3> a1 = polar(radius, a, z1);
    const Vector<T, 3> b1 = polar(radius, b, z1);
    tris.push_back({a0, b0, a1});
    tris.push_back({b0, b1, a1});
    const Vector<T, 3> c0{Center, Center, z0};
    const Vector<T, 3> c1{Center, Center, z1};
    tris.push_back({c0, b0, a0});
    tris.push_back({c1, a1, b1});
  }
}

// Annular solid between rIn and rOut: the inner side points towards the axis.
void addAnnularShell(std::vector<Tri>& tris, T rIn, T rOut, T z0, T z1) {
  for (int i = 0; i < ThetaSegments; ++i) {
    const T a = T(2) * PiConst * T(i) / T(ThetaSegments);
    const T b = T(2) * PiConst * T(i + 1) / T(ThetaSegments);
    const Vector<T, 3> i0a = polar(rIn, a, z0);
    const Vector<T, 3> i0b = polar(rIn, b, z0);
    const Vector<T, 3> i1a = polar(rIn, a, z1);
    const Vector<T, 3> i1b = polar(rIn, b, z1);
    const Vector<T, 3> o0a = polar(rOut, a, z0);
    const Vector<T, 3> o0b = polar(rOut, b, z0);
    const Vector<T, 3> o1a = polar(rOut, a, z1);
    const Vector<T, 3> o1b = polar(rOut, b, z1);
    tris.push_back({i0a, i1a, i0b});
    tris.push_back({i0b, i1a, i1b});
    tris.push_back({o0a, o0b, o1a});
    tris.push_back({o0b, o1b, o1a});
    tris.push_back({i0a, i0b, o0b});
    tris.push_back({i0a, o0b, o0a});
    tris.push_back({i1a, o1b, i1b});
    tris.push_back({i1a, o1a, o1b});
  }
}

void writeStl(const std::string& path, const std::vector<Tri>& tris) {
  std::ofstream out(path);
  out << std::setprecision(17);
  out << "solid mesh\n";
  for (const Tri& tri : tris) {
    const Vector<T, 3> n =
      CrossProduct3(tri[1] - tri[0], tri[2] - tri[0]).getnormalize();
    out << "facet normal " << n[0] << ' ' << n[1] << ' ' << n[2] << "\n";
    out << "  outer loop\n";
    for (const Vector<T, 3>& v : tri) {
      out << "    vertex " << v[0] << ' ' << v[1] << ' ' << v[2] << "\n";
    }
    out << "  endloop\n";
    out << "  endfacet\n";
  }
  out << "endsolid mesh\n";
}

T analytic(T r) {
  const T ri2 = InnerRadius * InnerRadius;
  const T ro2 = OuterRadius * OuterRadius;
  const T A = -Omega * ri2 / (ro2 - ri2);
  const T B = Omega * ri2 * ro2 / (ro2 - ri2);
  return A * r + B / r;
}

}  // namespace

int main() {
  const T h = T(1);
  const T z0 = -T(2);
  const T z1 = T(Nz) + T(2);
  std::vector<Tri> innerMesh;
  std::vector<Tri> outerMesh;
  addCylinder(innerMesh, InnerRadius, z0, z1);
  addAnnularShell(outerMesh, OuterRadius, ShellOuterRadius, z0, z1);
  writeStl("inner_cylinder.stl", innerMesh);
  writeStl("outer_shell.stl", outerMesh);

  const AABB<T, 3> domain(Vector<T, 3>{T(0), T(0), T(0)},
                          Vector<T, 3>{T(Nx), T(Ny), T(Nz)});
  flbgeo::VoxelGeometryBuilder<T> builder(Nx, Ny, Nz, domain, h, T(1.5));
  builder.addBox(domain, 1, 1, flbgeo::VoxelRole::Domain);
  StlReader<T> innerReader("inner_cylinder.stl", h, T(1), 2);
  StlReader<T> outerReader("outer_shell.stl", h, T(1), 2);
  builder.addStl(innerReader, 2, 2, flbgeo::VoxelRole::Solid, true);
  builder.addStl(outerReader, 3, 3, flbgeo::VoxelRole::Solid, true);
  const flbgeo::VoxelGeometry<T> record = builder.build();
  std::remove("inner_cylinder.stl");
  std::remove("outer_shell.stl");

  const std::size_t domainCells = record.countRole(flbgeo::VoxelRole::Domain);
  const std::size_t solidCells = record.countRole(flbgeo::VoxelRole::Solid);
  CHECK(domainCells > 1000);
  CHECK(solidCells > 0);
  CHECK(record.getBandSize() > 0);

  std::size_t innerSamples = 0;
  std::size_t outerSamples = 0;
  for (std::size_t id : record.getBandCells()) {
    if (!record.isDomain(id)) continue;
    if (record.getSurfaceSample(id).Entity == 2) ++innerSamples;
    if (record.getSurfaceSample(id).Entity == 3) ++outerSamples;
  }
  CHECK(innerSamples > 0);
  CHECK(outerSamples > 0);

  BlockGeometry3D<T> geo(Nx, Ny, Nz, 1, record.getDomain(),
                         record.getVoxelSize(), 1);
  BlockFieldManager<FLAG, T, 3> flagFM(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, flagFM, BulkFlag, VoidFlag);

  // Periodic z faces: the geometry is invariant along z, so the caps of both
  // meshes sit outside the domain and the profile has no end effects.
  const AABB<T, 3> bottom(Vector<T, 3>{T(0), T(0), T(0)},
                          Vector<T, 3>{T(Nx), T(Ny), T(1)});
  const AABB<T, 3> top(Vector<T, 3>{T(0), T(0), T(Nz - 1)},
                       Vector<T, 3>{T(Nx), T(Ny), T(Nz)});
  auto makePeriodic = [&](const AABB<T, 3>& box) {
    flagFM.forEach(box, [&](FLAG& field, std::size_t id) {
      if (util::isFlag(field.get(id), BulkFlag)) {
        field.SetField(id, std::uint8_t(BulkFlag | PeriodicFlag));
      }
    });
  };
  makePeriodic(bottom);
  makePeriodic(top);

  BaseConverter<T> conv(LatSet::cs2);
  conv.SimplifiedConverterFromRT(Nx, T(0.1), T(0.56));
  ValuePack initValues(conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(geo, initValues, conv);

  FixedPeriodicBoundaryManager<BlockLatticeManager<T, LatSet, FIELDS>,
                               BlockFieldManager<FLAG, T, 3>>
    NS_Peri("NS_Peri", NSLattice, flagFM, PeriodicFlag, VoidFlag);
  NS_Peri.Setup(bottom, NbrDirection::ZN, top, NbrDirection::ZP);

  std::unordered_map<std::uint16_t, flbbc::WallMotion<T>> motions;
  motions[2] = flbbc::WallMotion<T>{Vector<T, 3>{Center, Center, T(0)},
                                    Vector<T, 3>{T(0), T(0), T(1)}, Omega};

  flbbc::InterpolatedBBManager<decltype(NSLattice), decltype(flagFM)> ibb(
    NSLattice, flagFM, record, VoidFlag, PeriodicFlag, nullptr, &motions);
  CHECK(ibb.getLinkCount() > 0);

  std::size_t movingLinks = 0;
  std::size_t staticLinks = 0;
  for (const flbbc::SurfaceLink<T>& link : ibb.getLinks(0)) {
    if (link.WallVelocity.getnorm2() > T(0)) {
      ++movingLinks;
    } else {
      ++staticLinks;
    }
  }
  CHECK(movingLinks > 0);
  CHECK(staticLinks > 0);
  std::cout << "[feature/link_bc_curved] domain=" << domainCells
            << " solid=" << solidCells << " links=" << ibb.getLinkCount()
            << " inner_samples=" << innerSamples
            << " outer_samples=" << outerSamples
            << " moving_links=" << movingLinks
            << " static_links=" << staticLinks << std::endl;

  for (int step = 0; step < Steps; ++step) {
    NSLattice.ApplyCellDynamics<NSTask>(flagFM);
    NS_Peri.Apply();
    NSLattice.Stream();
    ibb.Apply();
  }
  NSLattice.ApplyCellDynamics<TaskSelectorRhoU>(flagFM);

  const auto& block = geo.getBlock(0);
  const int overlap = block.getOverlap();
  const int zm = Nz / 2;
  const auto& flag = flagFM.getBlockField(0).getField(0);
  const auto& vel = NSLattice.getBlockLat(0).template getField<VELOCITY<T, LatSet::d>>();

  T sqErr{};
  T sqErrInner{};
  T sqErrOuter{};
  std::size_t samples = 0;
  std::size_t samplesInner = 0;
  std::size_t samplesOuter = 0;
  T maxNum{};
  T maxRef{};
  const T midRadius = T(0.5) * (InnerRadius + OuterRadius);
  for (int y = 0; y < Ny; ++y) {
    for (int x = 0; x < Nx; ++x) {
      const std::size_t id = block.getIndex(
        Vector<int, 3>{x + overlap, y + overlap, zm + overlap});
      if (flag[id] == VoidFlag) continue;
      if (util::isFlag(flag[id], PeriodicFlag)) continue;
      const Vector<T, 3> pos = block.getVoxel(Vector<int, 3>{x + overlap,
                                                             y + overlap,
                                                             zm + overlap});
      const T xr = pos[0] - Center;
      const T yr = pos[1] - Center;
      const T r = std::sqrt(xr * xr + yr * yr);
      if (r <= InnerRadius || r >= OuterRadius) continue;
      const Vector<T, LatSet::d> u = vel.get(id);
      const T vtheta = (xr * u[1] - yr * u[0]) / r;
      const T ref = analytic(r);
      const T d = vtheta - ref;
      sqErr += d * d;
      if (r < midRadius) {
        sqErrInner += d * d;
        ++samplesInner;
      } else {
        sqErrOuter += d * d;
        ++samplesOuter;
      }
      maxNum = std::max(maxNum, std::abs(vtheta));
      maxRef = std::max(maxRef, std::abs(ref));
      ++samples;
    }
  }
  CHECK(samples > 100);
  CHECK(maxRef > T(0));
  const T l2 = std::sqrt(sqErr / T(samples)) / maxRef;
  const T l2Inner =
    std::sqrt(sqErrInner / T(std::max<std::size_t>(samplesInner, 1))) / maxRef;
  const T l2Outer =
    std::sqrt(sqErrOuter / T(std::max<std::size_t>(samplesOuter, 1))) / maxRef;
  std::cout << "[feature/link_bc_curved] samples=" << samples
            << " z_plane=" << zm << " max_u_theta=" << maxNum
            << " max_ref=" << maxRef << " profile_l2=" << l2
            << " inner_l2=" << l2Inner << " outer_l2=" << l2Outer << std::endl;

  CHECK(l2 < T(0.10));

  return flbtest::summary("feature/link_bc_curved");
}
