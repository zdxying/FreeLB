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

#pragma once

//freelb.h
//include haeder files for LBM
//
// NOTE: the legacy PopLattice path is not usable through this umbrella --
// data_struct/lattice.h(hh) (PopLattice/RhoLattice, method definitions),
// lbm/buoyancy.h, lbm/legacy_collision.h (legacy::BGK) and the PopLattice-path
// CA::ZhuStefanescu* in ca/zhu_stefanescu*_legacy.h(hh) are NOT wrapped;
// include them directly if you need them (examples/legacy/*).
// The active BlockLattice-path CA declarations (ca/zhu_stefanescu2d.h,
// wrapped by freelb.hh) no longer pull in data_struct/lattice.h.  The legacy
// path is pending removal together with examples/legacy.

// timer
#include "utils/timer.h"
// geometry
#include "geometry/geometry.h"
// boundary
#include "boundary/boundary.h"

// io
#include "io/vtkWriter.h"
#include "io/vtm_writer.h"
#include "io/ini_reader.h"
#include "io/vtu_writer.h"
// checkpoint / restart; implementation follows in freelb.hh
#include "io/checkpoint.h"

#include "lbm/lattice_set.h"
// block-field statistics (FieldStatistics; used by examples/dambreak*)
#include "data_struct/field_statistics.h"
// lbm dynamics
#include "lbm/lbm.h"

// ----------------- extended models -----------------
// (kept outside the core single-phase path but wrapped here, so examples only
// need freelb.h + freelb.hh; see README "Umbrella headers")
// free-surface model
#include "lbm/freeSurface.h"
// STL geometry to field conversion
#include "offLattice/marchingCube.h"
// step-wise 2D data writer
#include "io/plot_writer.h"
// field checksum (regression validation)
#include "utils/field_checksum.h"

