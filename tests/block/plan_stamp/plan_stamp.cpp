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

// plan_stamp.cpp
// Plan generation stamps: a communication plan must be judged against the
// geometry generation it was built at, not against whatever the manager built
// last.
//
// BlockLatticeManager used to keep one _planGeneration updated by every build,
// so building the cross-level plan after the geometry changed also marked the
// older full plan valid again.  Its cell indices then addressed the rebuilt
// geometry.  This target reproduces that sequence on one process:
//
//   1. build the full plan and record its generation,
//   2. bump the geometry generation,
//   3. build the Aver cross plan,
//   4. require the full plan to keep its old stamp, the cross plan to carry the
//      current one, and the full plan to differ from the current generation.
//
// The stamps are compared directly: the consumer path would abort (MPI_Abort)
// on a mismatch, and a test must not take the process down.
//
// One process is enough; the plan stamp does not depend on the rank count.

#include <cstdio>

#include "freelb.h"
#include "freelb.hh"
#include "../../test_util.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);
  if (mpi().getSize() != 1) {
    std::fprintf(stderr, "[plan_stamp] run this target on one rank\n");
    return 1;
  }

  AABB<T, 2> domain(Vector<T, 2>{}, Vector<T, 2>{T(32), T(32)});
  BlockGeometryHelper2D<T> helper(32, 32, domain, T(1), 16);
  helper.CreateBlocks();
  helper.LoadBalancing(1);
  BlockGeometry2D<T> geo(helper);

  BaseConverter<T> conv(LatSet::cs2);
  conv.ConvertFromRT(T(1), T(0.56), T(0.001), T(32), T(0.1), T(0.02));
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  ValuePack<T, Vector<T, LatSet::d>, T> initValues(conv.getLatRhoInit(),
                                                   Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> lattice(geo, initValues, conv);

  // Start above zero so a plan whose stamp was never assigned (the default 0)
  // cannot pass by accident.
  geo.bumpGeometryGeneration();
  const std::uint64_t built = geo.getGeometryGeneration();

  auto& fullPlan = lattice.planFor(flbcomm::CommMode::Full, 0);
  CHECK_EQ(fullPlan.Generation, built);

  geo.bumpGeometryGeneration();
  const std::uint64_t current = geo.getGeometryGeneration();
  CHECK(current != built);

  auto& crossPlan = lattice.crossPlanFor(flbcomm::CrossMode::Aver, 0);
  CHECK_EQ(crossPlan.Generation, current);
  CHECK_EQ(fullPlan.Generation, built);
  CHECK(fullPlan.Generation != current);

  std::printf("[plan_stamp] built=%llu current=%llu full=%llu cross=%llu\n",
              static_cast<unsigned long long>(built),
              static_cast<unsigned long long>(current),
              static_cast<unsigned long long>(fullPlan.Generation),
              static_cast<unsigned long long>(crossPlan.Generation));
  return flbtest::summary("block/plan_stamp");
}
