/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 */

// instantiation_smoke.cpp
//
// Compile-time smoke test for header self-containedness.
//
// Rationale: several lbm schemes are covered by NO example target.  If a header
// relies on transitive includes that later get cleaned up, every build stays
// green until someone actually uses the scheme.  This TU explicitly
// instantiates the uncovered schemes so include regressions become build
// errors here, instead of latent template-instantiation errors downstream.
//
// Covered here:
//   * collision::MRT_Feq_RhoU (MRT path -- zero coverage in examples/)
//   * equilibrium::SecondOrder, moment::rhoU, force::ForcePop (via MRT)
//   * CA::ZhuStefanescu3D full member instantiation (legacy/cazs3d covers its
//     used subset; this covers the rest)
//   * bounce-back / periodic boundary UpdateRho (moment::Rho)
// Covered by examples/ (dambreak2d/dambreak3d): free-surface schemes.
// The MPI dynamic-refinement path in data_struct/block_lattice.hh is forced by
// the `mpi` Makefile target (-DMPI_ENABLED); the default target is CPU-only.

#include "freelb.h"
#include "freelb.hh"

// Block CA is part of the active path.  freelb.hh already pulls this header
// in (ca/zhu_stefanescu2d.hh -> .h); it is included explicitly here so the
// test stays valid if the umbrella contract changes.
#include "ca/zhu_stefanescu2d.h"

#include "boundary/bounce_back_boundary.h"
#include "boundary/periodic_boundary.h"
#include "lbm/collisionMRT.h"

using T = FLOAT;
using LatSet = D3Q19<T>;
using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL = Cell<T, LatSet, FIELDS>;

// Explicit instantiations: each of these forces the full member bodies to be
// type-checked, including every equilibrium:: / moment:: / force:: reference.
template struct collision::MRT_Feq_RhoU<CELL>;
template struct collision::MRT_Feq_RhoU<CELL, true>;

// CA solidification, ACTIVE Block path (matches examples/cazsblock2d):
// full-member instantiation of the BlockLattice-based CA.
template class CA::BlockZhuStefanescu2D<T, CA::D2Q8<T>>;

// Non-Newtonian (power-law) and Smagorinsky LES schemes: zero consumers in
// examples/, so they would otherwise never be instantiated (and thus never
// type-checked).  Both need extra per-cell parameter fields, so build matching
// cell types and force the whole member bodies.
using NN_FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>,
                           MinShearRate<T>, MaxShearRate<T>, PL_k<T>, PL_m<T>,
                           OMEGA<T>>;
using NN_CELL = Cell<T, LatSet, NN_FIELDS>;
template struct PowerLaw_Omega<NN_CELL>;
template struct collision::Rheology_BGK<moment::rhoU<NN_CELL>,
                                        equilibrium::SecondOrder<NN_CELL>,
                                        PowerLaw_Omega<NN_CELL>>;

using SMAG_FIELDS =
    TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>, SMAGORINSKY<T>>;
using SMAG_CELL = Cell<T, LatSet, SMAG_FIELDS>;
template struct collision::SmagorinskyBGK<moment::rhoU<SMAG_CELL>,
                                          equilibrium::SecondOrder<SMAG_CELL>>;

#ifdef MPI_ENABLED
// Dynamic block refinement/coarsening driver (block_lattice.hh).  Its whole
// method set is forced here because no example drives it (the legacy AMR
// examples were removed).  Build this TU with -DMPI_ENABLED (see Makefile
// `mpi` target) to exercise the MPI-guarded code paths inside block_lattice.hh.
// The helper is 2D-only (it uses BlockGeometry2D), so instantiate it with a
// 2D lattice set.
using LatSet2D = D2Q9<T>;
using DYN_FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet2D::d>, POP<T, LatSet2D::q>>;
template class DynamicBlockLatticeHelper2D<T, LatSet2D, DYN_FIELDS>;
#endif

int main() { return 0; }
