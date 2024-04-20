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

// basic_boundary.h

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>

#include "data_struct/block_lattice.h"
#include "lbm/equilibrium.h"
#include "utils/alias.h"

// fixed boundary cell structure
struct FixedBdCell {
  // outflow directions
  std::vector<int> outflows;
  // cell id
  std::size_t Id;

  FixedBdCell(std::size_t id, int q) : Id(id) { outflows.reserve(q); }
};

class AbstractBoundary {
 public:
  virtual void Apply() = 0;
  virtual void getinfo() {}
  virtual void UpdateRho() {}
  virtual void UpdateU() {}
  // virtual void UpdateU_F() = 0;
};

// flagType should be std::uint8_t or enum of std::uint8_t
template <typename T, typename LatSet, typename flagType>
class FixedBoundary : public AbstractBoundary {
 protected:
  // boundary cell
  std::vector<FixedBdCell> BdCells;
  // reference to lattice
  BasicLattice<T, LatSet> &Lat;
  // reference to geometry
  Geometry<T, LatSet::d> &Geo;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  // geometry flag
  GenericArray<flagType> &Field;

 public:
  FixedBoundary(BasicLattice<T, LatSet> &lat, std::uint8_t cellflag,
                std::uint8_t voidflag);
  FixedBoundary(BasicLattice<T, LatSet> &lat, GenericArray<flagType> &f,
                std::uint8_t cellflag, std::uint8_t voidflag);
  // get boundary cell flag
  std::uint8_t getBdCellFlag() const { return BdCellFlag; }
  // get void cell flag
  std::uint8_t getVoidFlag() const { return voidFlag; }
  // add to boundary cells
  void addtoBd(std::size_t id);
  // setup boundary cells
  void Setup();
};

template <typename T, typename LatSet, typename flagType>
class MovingBoundary : public AbstractBoundary {
 protected:
  // boundary cells
  std::vector<std::size_t> &Ids;
  // reference to lattice
  BasicLattice<T, LatSet> &Lat;
  // reference to geometry
  Geometry<T, LatSet::d> &Geo;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  // geometry flag
  GenericArray<flagType> &Field;

 public:
  MovingBoundary(BasicLattice<T, LatSet> &lat, std::vector<std::size_t> &ids,
                 std::uint8_t voidflag, std::uint8_t cellflag);
  MovingBoundary(BasicLattice<T, LatSet> &lat, std::vector<std::size_t> &ids,
                 GenericArray<flagType> &f, std::uint8_t voidflag, std::uint8_t cellflag);

  // get boundary cell flag
  std::uint8_t getBdCellFlag() const { return BdCellFlag; }
  // get void cell flag
  std::uint8_t getVoidFlag() const { return voidFlag; }
  // get std::vector<std::size_t> &Ids;
  std::vector<std::size_t> &getIds() { return Ids; }
  // update boundary cells: std::vector<std::size_t> &Ids;
  void UpdateBdCells();
};

// --------------------------------------------------------------------------------------
// -----------------------------------BlockBoundary--------------------------------------
// --------------------------------------------------------------------------------------

class AbstractBlockBoundary {
 public:
  virtual void Apply(std::int64_t count) = 0;
  virtual void UpdateRho(std::int64_t count) {}
  virtual void UpdateU(std::int64_t count) {}
};

template <typename T, typename LatSet, typename flagType>
class BlockFixedBoundary {
 protected:
  // boundary cell
  std::vector<FixedBdCell> BdCells;
  // reference to lattice
  BlockLattice<T, LatSet> &Lat;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  // geometry flag
  const GenericArray<flagType> &Field;

 public:
  BlockFixedBoundary(BlockLattice<T, LatSet> &lat, const GenericArray<flagType> &f,
                     std::uint8_t cellflag, std::uint8_t voidflag);
  // get boundary cell flag
  std::uint8_t getBdCellFlag() const { return BdCellFlag; }
  // get void cell flag
  std::uint8_t getVoidFlag() const { return voidFlag; }
  BlockLattice<T, LatSet> &getLat() { return Lat; }
  // add to boundary cells: std::vector<FixedBdCell> BdCells
  void addtoBd(std::size_t id);
  // setup boundary cells
  void Setup();
};

template <typename T, typename LatSet, typename flagType>
class BlockMovingBoundary {
 protected:
  // boundary cells
  std::vector<std::size_t> &Ids;
  // reference to lattice
  BlockLattice<T, LatSet> &Lat;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  // geometry flag
  GenericArray<flagType> &Field;

 public:
  BlockMovingBoundary(BlockLattice<T, LatSet> &lat, std::vector<std::size_t> &ids,
                      GenericArray<flagType> &f, std::uint8_t voidflag,
                      std::uint8_t cellflag);
  // get boundary cell flag
  std::uint8_t getBdCellFlag() const { return BdCellFlag; }
  // get void cell flag
  std::uint8_t getVoidFlag() const { return voidFlag; }
  // get boundary cells std::vector<std::size_t> &Ids;
  std::vector<std::size_t> &getIds() { return Ids; }
  BlockLattice<T, LatSet> &getLat() { return Lat; }
  // update boundary cells
  void UpdateBdCells();
};

class BoundaryManager {
 private:
  std::vector<AbstractBoundary *> _Boundaries;

 public:
  BoundaryManager(std::vector<AbstractBoundary *> boundaries) : _Boundaries(boundaries) {
    printinfo();
  }
  template <typename... Args>
  BoundaryManager(Args... args) : _Boundaries{args...} {
    printinfo();
  }

  void Apply() {
    for (AbstractBoundary *boundary : _Boundaries) boundary->Apply();
  }
  void printinfo() {
    std::cout << "[Boundary Statistics]: "
              << "\n"
              << "Boundary Type  |  Number of Boundary Cells" << std::endl;
    for (AbstractBoundary *boundary : _Boundaries) boundary->getinfo();
  }
  void UpdateRho() {
    for (AbstractBoundary *boundary : _Boundaries) boundary->UpdateRho();
  }
  void UpdateU() {
    for (AbstractBoundary *boundary : _Boundaries) boundary->UpdateU();
  }
};

class BlockBoundaryManager {
 private:
  std::vector<AbstractBlockBoundary *> _Boundaries;

 public:
  BlockBoundaryManager(std::vector<AbstractBlockBoundary *> boundaries)
      : _Boundaries(boundaries) {}

  template <typename... Args>
  BlockBoundaryManager(Args... args) : _Boundaries{args...} {}

  void Apply(std::int64_t count) {
    for (AbstractBlockBoundary *boundary : _Boundaries) boundary->Apply(count);
  }
  void UpdateRho(std::int64_t count) {
    for (AbstractBlockBoundary *boundary : _Boundaries) boundary->UpdateRho(count);
  }
  void UpdateU(std::int64_t count) {
    for (AbstractBlockBoundary *boundary : _Boundaries) boundary->UpdateU(count);
  }
};
