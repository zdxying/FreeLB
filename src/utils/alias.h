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
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

// alias.h

#pragma once

#include <array>
#include <cstdint>
#include <type_traits>

template <unsigned int D>
struct BasicInterp;
using BasicInterp2D = BasicInterp<4>;
using BasicInterp3D = BasicInterp<8>;
template <unsigned int D>
using InterpStruct = std::conditional_t<D == 2, BasicInterp2D, BasicInterp3D>;


using InterpSource2D = std::array<std::size_t, 4>;
using InterpSource3D = std::array<std::size_t, 8>;
// std::array<int, 4/8>
template <unsigned int D>
using InterpSource = std::conditional_t<D == 2, InterpSource2D, InterpSource3D>;

template <typename T>
using InterpWeight2D = std::array<T, 4>;
template <typename T>
using InterpWeight3D = std::array<T, 8>;
// std::array<T, 4/8>
template <typename T, unsigned int D>
using InterpWeight = std::conditional_t<D == 2, InterpWeight2D<T>, InterpWeight3D<T>>;


// ---------geometry alias-----------
template <typename T, unsigned int D>
class AABB;

template <typename T>
class Geometry2D;
template <typename T>
class Geometry3D;
template <typename T, unsigned int D>
using Geometry = std::conditional_t<D == 2, Geometry2D<T>, Geometry3D<T>>;

template <typename T>
class Block2D;
template <typename T>
class Block3D;
template <typename T, unsigned int D>
using Block = std::conditional_t<D == 2, Block2D<T>, Block3D<T>>;

template <typename T>
class BlockGeometry2D;
template <typename T>
class BlockGeometry3D;
template <typename T, unsigned int D>
using BlockGeometry = std::conditional_t<D == 2, BlockGeometry2D<T>, BlockGeometry3D<T>>;

template <typename T>
class VoxelGeometry2D;
template <typename T>
class VoxelGeometry3D;
template <typename T, unsigned int D>
using VoxelGeometry = std::conditional_t<D == 2, VoxelGeometry2D<T>, VoxelGeometry3D<T>>;

template <typename T>
class RefinedGeometry2D;
template <typename T>
class RefinedGeometry3D;
template <typename T, unsigned int D>
using RefinedGeometry = std::conditional_t<D == 2, RefinedGeometry2D<T>, RefinedGeometry3D<T>>;

// ---------field alias-----------
template <typename ArrayType, unsigned int D>
class GenericField;

template <typename T>
class GenericArray;

template <typename T>
class CyclicArray;

template <typename FieldType, typename T>
class BlockFieldStruct;

using FlagArray = GenericArray<std::uint8_t>;
using FlagField = GenericField<FlagArray, 1>;

template <typename T>
using ScalerField = GenericField<GenericArray<T>, 1>;

// array of structure version of vector field
// access: getField()[index][ith component]
template <typename T, unsigned int D>
using VectorFieldAOS = GenericField<GenericArray<Vector<T, D>>, 1>;

// structure of array version of vector field
// access: getField(ith component)[index]
template <typename T, unsigned int D>
using VectorFieldSoA = GenericField<GenericArray<T>, D>;

template <typename T, unsigned int q>
using PopulationField = GenericField<CyclicArray<T>, q>;

template <typename T, unsigned int D>
using BlockVectFieldAOS = BlockFieldStruct<VectorFieldAOS<T, D>, Vector<T, D>>;

template <typename T>
using BlockScalerField = BlockFieldStruct<ScalerField<T>, T>;

template <template <typename> class ArrayType, typename T, unsigned int D>
using BlockFStruct = BlockFieldStruct<GenericField<ArrayType<T>, D>, T>;

namespace CA {
template <typename T, typename LatSet>
class BlockZhuStefanescu2D;
template <typename T, typename LatSet>
class BlockZhuStefanescu3D;
template <typename T, typename LatSet>
using BlockZhuStefanescu =
  std::conditional_t<LatSet::d == 2, BlockZhuStefanescu2D<T, LatSet>, BlockZhuStefanescu3D<T, LatSet>>;
}  // namespace CA
