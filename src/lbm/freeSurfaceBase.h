/* This file is part of FreeLB, modified from openLB and FluidX3D with the following copyright notice:
 *
 * // start of the original OpenLB's copyright notice
 * 
 * This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 * 
 * // end of the original OpenLB's copyright notice
 * 
 * FluidX3D: https://github.com/ProjectPhysX/FluidX3Ds
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

// base of free surface model
#pragma once

/* This file is part of FreeLB, modified from openLB and FluidX3D with the following copyright notice:
 *
 * // start of the original OpenLB's copyright notice
 * 
 * This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 * 
 * // end of the original OpenLB's copyright notice
 * 
 * FluidX3D: https://github.com/ProjectPhysX/FluidX3Ds
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

// free surface model
#pragma once

#include "utils/alias.h"


// namespace olbfs: openlb's implementation of free surface model

namespace olbfs {

// free surface cell type
enum FSType : std::uint8_t {
  Void = 1,
  Wall = 2,
  Gas = 4,
  Interface = 8,
  Fluid = 16
};

// free surface transition flag
enum FSFlag : std::uint8_t {
  None = 0,
  To_Fluid = 1,
  To_Gas = 2,
  To_Interface = 4
};

// define unique olbfs Field
struct STATEBase : public FieldBase<1> {};
struct FLAGBase : public FieldBase<1> {};
struct MASSBase : public FieldBase<1> {};
struct VOLUMEFRACBase : public FieldBase<1> {};
struct MASSEXBase : public FieldBase<1> {};
struct PREVIOUS_VELOCITYBase : public FieldBase<1> {};

// free surface state, init with Solid
using STATE = GenericField<GenericArray<FSType>, STATEBase>;
// free surface transition flag
using FLAG = GenericField<GenericArray<FSFlag>, FLAGBase>;
// mass = rho * volumefraction
template <typename T>
using MASS = GenericField<GenericArray<T>, MASSBase>;
// fill level/ volume fraction in VOF
template <typename T>
using VOLUMEFRAC = GenericField<GenericArray<T>, VOLUMEFRACBase>;
// a simple massex scalar filed is enough for excess mass exchange
// more efficient than using a vector field of size q*N in openlb
template <typename T, unsigned int q>
using MASSEX = GenericField<GenericArray<Vector<T, q>>, MASSEXBase>;
// previous velocity in openlb
template <typename T, unsigned int D>
using PREVIOUS_VELOCITY = GenericField<GenericArray<Vector<T, D>>, PREVIOUS_VELOCITYBase>;


// define olbfs parameters as single data stored in Data
struct Lonely_ThBase : public FieldBase<1> {};
struct VOF_Trans_ThBase : public FieldBase<1> {};
struct Surface_Tension_EnabledBase : public FieldBase<1> {};
struct Surface_Tension_ParameterBase : public FieldBase<1> {};

// lonely threshold in mass transfer
template <typename T>
using Lonely_Th = Data<T, Lonely_ThBase>;
// vof transition threshold
template <typename T>
using VOF_Trans_Th = Data<T, VOF_Trans_ThBase>;
// surface tension enabled
using Surface_Tension_Enabled = Data<bool, Surface_Tension_EnabledBase>;
// surface tension parameter
template <typename T>
using Surface_Tension_Parameter = Data<T, Surface_Tension_ParameterBase>;


template <typename T, typename LatSet>
using FSFIELDS = TypePack<STATE, FLAG, MASS<T>, VOLUMEFRAC<T>, MASSEX<T, LatSet::q>, PREVIOUS_VELOCITY<T, LatSet::d>>;

template <typename T>
using FSPARAMS = TypePack<Lonely_Th<T>, VOF_Trans_Th<T>, Surface_Tension_Enabled,
                          Surface_Tension_Parameter<T>>;

}  // namespace olbfs


// FluidX3D's implementation of free surface model

namespace fx3dfs {

enum FSType : std::uint8_t {
  Solid = 1,
  Wall = 2,
  Gas = 4,
  Interface = 8,
  Fluid = 16,
  To_Fluid = 32,
  To_Gas = 64,
  To_Interface = 128
};


struct STATEBase : public FieldBase<1> {};
struct MASSBase : public FieldBase<1> {};
struct VOLUMEFRACBase : public FieldBase<1> {};
struct MASSEXBase : public FieldBase<1> {};

// free surface state, init with Solid
using STATE = GenericField<GenericArray<FSType>, STATEBase>;
// mass = rho * volumefraction
template <typename T>
using MASS = GenericField<GenericArray<T>, MASSBase>;
// fill level/ volume fraction in VOF
template <typename T>
using VOLUMEFRAC = GenericField<GenericArray<T>, VOLUMEFRACBase>;
// Excess mass, Lehmann's implementation
template <typename T>
using MASSEX = GenericField<GenericArray<T>, MASSEXBase>;


struct VOF_Trans_ThBase : public FieldBase<1> {};
struct Surface_Tension_EnabledBase : public FieldBase<1> {};
struct Surface_Tension_ParameterBase : public FieldBase<1> {};

// vof transition threshold
template <typename T>
using VOF_Trans_Th = Data<T, VOF_Trans_ThBase>;
// surface tension enabled
using Surface_Tension_Enabled = Data<bool, Surface_Tension_EnabledBase>;
// surface tension parameter
template <typename T>
using Surface_Tension_Parameter = Data<T, Surface_Tension_ParameterBase>;


template <typename T, typename LatSet>
using FSFIELDS = TypePack<STATE, MASS<T>, VOLUMEFRAC<T>, MASSEX<T>>;

template <typename T>
using FSPARAMS = TypePack<VOF_Trans_Th<T>, Surface_Tension_Enabled, Surface_Tension_Parameter<T>>;

}  // namespace fx3dfs
