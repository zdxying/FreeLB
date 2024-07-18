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

// alias.h

#pragma once

// std::array
#include <array>
// size_t, uint8_t
#include <cstdint>

#include "head.h"

template <unsigned int D>
struct BasicInterp;
using BasicInterp2D = BasicInterp<4>;
using BasicInterp3D = BasicInterp<8>;
template <unsigned int D>
using InterpStruct = std::conditional_t<D == 2, BasicInterp2D, BasicInterp3D>;


using InterpSource2D = std::array<std::size_t, 4>;
using InterpSource3D = std::array<std::size_t, 8>;
// std::array<std::size_t, 4/8>
template <unsigned int D>
using InterpSource = std::conditional_t<D == 2, InterpSource2D, InterpSource3D>;

template <typename T>
using InterpWeight2D = std::array<T, 4>;
template <typename T>
using InterpWeight3D = std::array<T, 8>;
// std::array<T, 4/8>
template <typename T, unsigned int D>
using InterpWeight = std::conditional_t<D == 2, InterpWeight2D<T>, InterpWeight3D<T>>;

namespace cudev {
#ifdef __CUDACC__

using InterpSource2D = thrust::tuple<std::size_t, std::size_t, std::size_t, std::size_t>;
using InterpSource3D = thrust::tuple<std::size_t, std::size_t, std::size_t, std::size_t,
                                     std::size_t, std::size_t, std::size_t, std::size_t>;
template <unsigned int D>
using InterpSource = std::conditional_t<D == 2, InterpSource2D, InterpSource3D>;

template <typename T>
using InterpWeight2D = thrust::tuple<T, T, T, T>;
template <typename T>
using InterpWeight3D = thrust::tuple<T, T, T, T, T, T, T, T>;
template <typename T, unsigned int D>
using InterpWeight = std::conditional_t<D == 2, InterpWeight2D<T>, InterpWeight3D<T>>;

template <typename T, typename Base>
class Data;

template <typename T, typename Base>
class Array;

template <typename ArrayType, unsigned int D>
class GenericArrayField;

template <typename ArrayType, typename Base>
class GenericField;

template <typename T>
class GenericArray;

template <typename T>
class CyclicArray;

template <typename T>
class StreamArray;

#endif
}  // namespace cudev

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
class BlockGeometryHelper2D;
template <typename T>
class BlockGeometryHelper3D;
template <typename T, unsigned int D>
using BlockGeometryHelper =
  std::conditional_t<D == 2, BlockGeometryHelper2D<T>, BlockGeometryHelper3D<T>>;


// ---------field alias-----------
// field base class
template <unsigned int D>
struct FieldBase {
  static constexpr unsigned int array_dim = D;
};

template <typename T, typename Base>
class Data;

template <typename T, typename Base>
class Array;

template <typename ArrayType, unsigned int D>
class GenericArrayField;

template <typename ArrayType, typename Base>
class GenericField;

template <typename T>
class GenericArray;

template <typename T>
class CyclicArray;

template <typename T>
class StreamArray;


template <typename T, unsigned int D>
class Vector;

template <typename T>
using ScalarField = GenericArrayField<GenericArray<T>, 1>;

using FlagField = ScalarField<std::uint8_t>;

// array of structure version of vector field
// access: getField()[index][ith component]
template <typename T, unsigned int D>
using VectorFieldAOS = GenericArrayField<GenericArray<Vector<T, D>>, 1>;

// structure of array version of vector field
// access: getField(ith component)[index]
template <typename T, unsigned int D>
using VectorFieldSoA = GenericArrayField<GenericArray<T>, D>;

template <typename T, unsigned int q>
using PopulationField = GenericArrayField<StreamArray<T>, q>;
// using PopulationField = GenericArrayField<CyclicArray<T>, q>;


// specific field name for access by Cell interface, not alias
struct RHOBase : public FieldBase<1> {};
struct TEMPBase : public FieldBase<1> {};
struct CONCBase : public FieldBase<1> {};
struct VELOCITYBase : public FieldBase<1> {};
struct FLAGBase : public FieldBase<1> {};
struct FORCEBase : public FieldBase<1> {};
struct SCALARFORCEBase : public FieldBase<1> {};
struct CONSTFORCEBase : public FieldBase<1> {};
struct SCALARCONSTFORCEBase : public FieldBase<1> {};
struct CONSTRHOBase : public FieldBase<1> {};
struct CONSTUBase : public FieldBase<1> {};

struct RHOINITBase : public FieldBase<1> {};
struct TEMPINITBase : public FieldBase<1> {};
struct CONCINITBase : public FieldBase<1> {};
struct GBETABase : public FieldBase<1> {};
template <unsigned int q>
struct POPBase : public FieldBase<q> {};


template <typename T>
using RHO = GenericField<GenericArray<T>, RHOBase>;
template <typename T>
using TEMP = GenericField<GenericArray<T>, TEMPBase>;
template <typename T>
using CONC = GenericField<GenericArray<T>, CONCBase>;

template <typename T, unsigned int D>
using VELOCITY = GenericField<GenericArray<Vector<T, D>>, VELOCITYBase>;

using FLAG = GenericField<GenericArray<std::uint8_t>, FLAGBase>;

template <typename T, unsigned int D>
using FORCE = GenericField<GenericArray<Vector<T, D>>, FORCEBase>;

template <typename T>
using SCALARFORCE = GenericField<GenericArray<T>, SCALARFORCEBase>;

template <typename T, unsigned int D>
using CONSTFORCE = Data<Vector<T, D>, CONSTFORCEBase>;

template <typename T>
using SCALARCONSTFORCE = Data<T, SCALARCONSTFORCEBase>;

template <typename T, unsigned int q>
using POP = GenericField<StreamArray<T>, POPBase<q>>;
// using POP = GenericField<CyclicArray<T>, POPBase<q>>;

template <typename T>
using RHOINIT = Data<T, RHOINITBase>;
template <typename T>
using TEMPINIT = Data<T, TEMPINITBase>;
template <typename T>
using CONCINIT = Data<T, CONCINITBase>;

template <typename T>
using GBETA = Data<T, GBETABase>;

template <typename T>
using CONSTRHO = Data<T, CONSTRHOBase>;

template <typename T, unsigned int D>
using CONSTU = Data<Vector<T, D>, CONSTUBase>;

// ---------block field alias-----------
template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockField;

template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockFieldManager;


namespace CA {
template <typename T, typename LatSet>
class BlockZhuStefanescu2D;
template <typename T, typename LatSet>
class BlockZhuStefanescu3D;
template <typename T, typename LatSet>
using BlockZhuStefanescu =
  std::conditional_t<LatSet::d == 2, BlockZhuStefanescu2D<T, LatSet>,
                     BlockZhuStefanescu3D<T, LatSet>>;
}  // namespace CA
