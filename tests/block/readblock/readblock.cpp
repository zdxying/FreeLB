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
using LatSet = D3Q19<T>;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/

int main() {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  // constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  // constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
  // constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  // ------------------ define geometry ------------------
  BlockReader<T,2> blockreader2d("block2d");
  BlockReader<T,3> blockreader3d("block3d");

  BlockGeometry2D<T> Geo2d(blockreader2d);
  BlockGeometry3D<T> Geo3d(blockreader3d);

  BlockFieldManager<FlagField, T, 2> FlagFM2d(Geo2d, VoidFlag);
  BlockFieldManager<FlagField, T, 3> FlagFM3d(Geo3d, VoidFlag);

  // FlagFM.forEach(cavity, [&](FlagField& field, std::size_t id) { field.SetField(id, AABBFlag); });
  // FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);

  vtmwriter::ScalarWriter Writer2d("flag", FlagFM2d);
  vtmwriter::vtmWriter<T, 2> GeoWriter2d("Geo2", Geo2d);
  GeoWriter2d.addWriterSet(Writer2d);
  GeoWriter2d.WriteBinary();

  vtmwriter::ScalarWriter Writer3d("flag", FlagFM3d);
  vtmwriter::vtmWriter<T, 3> GeoWriter3d("Geo3", Geo3d);
  GeoWriter3d.addWriterSet(Writer3d);
  GeoWriter3d.WriteBinary();
}