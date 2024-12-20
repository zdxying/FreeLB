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

#include "freelb.h"
#include "freelb.hh"
#include "offLattice/marchingCube.hh"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Ni;
int Nj;
int Nk;
T Cell_Len;
int Thread_Num;
std::string work_dir;

template <typename T>
using ScalarBFM = BlockFieldManager<RHO<T>, T, 3>;

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t AABBFlag = std::uint8_t(2);
constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
constexpr std::uint8_t FluidFlag = std::uint8_t(8);

void readParam() {
  iniReader param_reader("mcpa.ini");
  work_dir = param_reader.getValue<std::string>("MC", "workdir_");
  Thread_Num = param_reader.getValue<int>("MC", "thread_num");
  Ni = param_reader.getValue<int>("MC", "Ni");
  Nj = param_reader.getValue<int>("MC", "Nj");
  Nk = param_reader.getValue<int>("MC", "Nk");
  Cell_Len = param_reader.getValue<T>("MC", "Cell_Len");
}

void setField(ScalarBFM<T>& BFM, const BlockFieldManager<FLAG, T, 3>& Flag) {
  // set interface cells
  BFM.forEach(Flag, [&](auto& field, auto& flagf, std::size_t id) {
    if (util::isFlag(flagf.get(id), FluidFlag)) {
      const auto& block = field.getBlock();
      for (int i = 1; i < LatSet::q; ++i) {
        std::size_t idn = id + latset::c<LatSet>(i) * block.getProjection();
        Vector<T, LatSet::d> loc_t = block.getLoc_t(idn);
        if (block.IsInside(loc_t)) {
          if (util::isFlag(flagf.get(idn), AABBFlag)) {
            field.get(idn) = T{0.5};
          }
        }
      }
    }
  });
}

int main() {
  readParam();

  // define geometry
  AABB<T, 3> cavity(Vector<T, 3>{},
                    Vector<T, 3>(T(Ni * Cell_Len), T(Nj * Cell_Len), T(Nk * Cell_Len)));
  AABB<T, 3> fluid(Vector<T, 3>{Cell_Len},
                   Vector<T, 3>(T(int(Ni / 2) * Cell_Len), T(int(Nj / 2) * Cell_Len),
                                T(int(Nk / 2) * Cell_Len)));
  BlockGeometry3D<T> Geo(Ni, Nj, Nk, Thread_Num, cavity, Cell_Len);

  BlockFieldManager<FLAG, T, 3> flag(Geo, VoidFlag);
  flag.forEach(cavity,
               [&](FLAG& field, std::size_t id) { field.SetField(id, FluidFlag); });
  flag.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  flag.forEach(fluid,
               [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });

  // scalar field
  ScalarBFM<T> scalar(Geo, T{});
  // scalar.forEach(fluid, [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
	scalar.forEach(flag, [&](auto& field, auto& flagf, std::size_t id) {
		if (util::isFlag(flagf.get(id), FluidFlag)) field.get(id) = T{1};
	});
  scalar.template SetupBoundary<LatSet>(cavity, T{});
  setField(scalar, flag);

  FieldStatistics stat(scalar);
  stat.printValueStatistics();

  // vector field
  BlockFieldManager<VELOCITY<T, 3>, T, 3> velo(Geo, Vector<T, 3>{T{1}, T{2}, T{3}});

  vtmwriter::ScalarWriter FlagWriter("flag", flag);
  vtmwriter::ScalarWriter ScalarWriter("scalar", scalar);
  vtmwriter::vtmWriter<T, 3> Writer("mc", Geo);
  Writer.addWriterSet(FlagWriter, ScalarWriter);
  Writer.WriteBinary();

	// ISO surface stl
  offlat::MarchingCubeSurface<T, RHO<T>> mc(scalar, T{0.5});
  offlat::TriangleSet<T> triangles;
  mc.generateIsoSurface(triangles);

  triangles.writeBinarySTL("marchingCube");

  vtuwriter::ScalarWriter vtuScalarWriter("scalar", scalar, triangles);
  vtuwriter::VectorWriter vtuVectorWriter("velocity", velo, triangles);
  vtuwriter::vtuManager<T> vtuManager("vtu", triangles);
  vtuManager.addWriter(vtuScalarWriter, vtuVectorWriter);
  vtuManager.Write();

  return 0;
}