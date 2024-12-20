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

//freelb.h
//include haeder files for LBM

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

#include "lbm/lattice_set.h"
// lbm dynamics
#include "lbm/lbm.h"

