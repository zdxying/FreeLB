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


using T = FLOAT;
// using LatSet = D2Q9<T>;
using LatSet = D3Q19<T>;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
int Nk;
T Cell_Len;
int BlockNum;

void readParam() {
  iniReader param_reader("divblock.ini");
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Nk = param_reader.getValue<int>("Mesh", "Nk");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockNum = param_reader.getValue<int>("Mesh", "BlockNum");
}

int main() {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
  constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  // Printer::Print_BigBanner(std::string("Initializing..."));  

  readParam();

  // ------------------ define geometry ------------------
  // AABB<T, 2> cavity(Vector<T, 2>{}, Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 3> cavity(Vector<T, 3>{},
                      Vector<T, 3>(T(Ni * Cell_Len), T(Nj * Cell_Len), T(Nk * Cell_Len)));
  // test geocomm
  Circle<T> circle(T(Ni * Cell_Len / 8),
                   Vector<T, 2>(T(Ni * Cell_Len / 2), T(Nj * Cell_Len / 2)));

  // BlockGeometry2D<T> Geo(Ni, Nj, BlockNum, cavity, Cell_Len);
  BlockGeometryHelper3D<T> GeoHelper(Ni, Nj, Nk, cavity, Cell_Len, Ni/5);
  GeoHelper.CreateBlocks(BlockNum);
  GeoHelper.AdaptiveOptimization(BlockNum);
  GeoHelper.LoadBalancing();

  // BlockGeometry3D<T> Geo(Ni, Nj, Nk, BlockNum, cavity, Cell_Len);
  BlockGeometry3D<T> Geo(GeoHelper);
  T dev = ComputeStdDev(Geo.getBasicBlocks());
  std::cout << "Standard deviation: " << dev << std::endl;

  BlockFieldManager<FlagField, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity,
                 [&](FlagField& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);

  BlockFieldManager<FlagField, T, LatSet::d> FlagFM2(Geo, VoidFlag);

  vtmwriter::ScalarWriter GeoFlagWriter("flag", FlagFM);
  vtmwriter::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo);
  GeoWriter.addWriterSet(GeoFlagWriter);

  GeoWriter.WriteBinary();
}