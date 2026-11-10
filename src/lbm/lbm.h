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

// lbm.h

// base schemes: loop implementations, always compiled
#include "lbm/moment.h"
#include "lbm/equilibrium.h"
#include "lbm/force.h"
#include "lbm/collision.h"
#include "lbm/collisionLES.h"
#include "lbm/collisionMRT.h"
#include "lbm/nonNewtonian.h"

// _FLB_UR_CSE is defined by make.mk iff -D_UNROLLFOR is given AND the csegen
// tool exists; the includes below then resolve to the cse-generated output in
// generated/lbm/.  In every other case (macro absent, or csegen missing) the
// base schemes above are used as-is: the hand-written copies under
// src/ur_fallback/lbm/ are a frozen comparison reference and are NEVER on the
// include path.
#ifdef _FLB_UR_CSE
// cse-generated unrolled specializations
#include "lbm/moment.ur.h"
#include "lbm/equilibrium.ur.h"
#include "lbm/force.ur.h"
// hand-written unrolled specializations without a cse counterpart
// (bounce-back collision / MRT)
#include "lbm/collision_unroll.h"
#include "lbm/collisionMRT_unroll.h"
#endif
