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

// plan_runner.h

#pragma once

#ifdef MPI_ENABLED

#include <cstddef>
#include <vector>

// Thread_Num
#include "head.h"
#include "parallel/comm_plan.h"
#include "parallel/mpi_manager.h"

namespace flbcomm {

// Parallel packing/unpacking pays off only when there is real work per group.
// Small plans (one peer, or a few hundred cells) are faster serially because
// the OpenMP region and the wait-all cost more than they save.
inline constexpr int kParallelMinGroups = 2;
inline constexpr std::size_t kParallelMinElements = 2048;

// Post every receive of a plan.  Receives go first so the peers can make
// progress while this rank is still packing its sends.
template <typename Plan>
void postRecvs(Plan& plan, std::vector<MPI_Request>& recvRequests) {
  recvRequests.resize(plan.Recv.size());
  for (std::size_t g = 0; g < plan.Recv.size(); ++g) {
    auto& group = plan.Recv[g];
    mpi().iRecv(group.Buffer.data(), static_cast<int>(group.Buffer.size()),
                group.TargetRank, &recvRequests[g], group.Tag);
  }
}

// Pack every send group, then post the sends.  Packing is the expensive part
// and touches only disjoint group buffers, so it runs in parallel; MPI itself
// stays on the calling thread (MPI_THREAD_FUNNELED).
template <typename Plan, typename PackFn>
void postSends(Plan& plan, std::vector<MPI_Request>& sendRequests, PackFn&& pack) {
  sendRequests.resize(plan.Send.size());
  const int groups = static_cast<int>(plan.Send.size());
  std::size_t total = 0;
  for (const auto& group : plan.Send) total += group.Buffer.size();
  auto packGroup = [&](int g) {
    auto& group = plan.Send[g];
    for (const CommRef& ref : group.Refs) {
      pack(ref, group.Buffer, ref.Offset);
    }
  };
  if (groups >= kParallelMinGroups && total >= kParallelMinElements) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
    for (int g = 0; g < groups; ++g) packGroup(g);
  } else {
    for (int g = 0; g < groups; ++g) packGroup(g);
  }
  for (std::size_t g = 0; g < plan.Send.size(); ++g) {
    auto& group = plan.Send[g];
    mpi().iSend(group.Buffer.data(), static_cast<int>(group.Buffer.size()),
                group.TargetRank, &sendRequests[g], group.Tag);
  }
}

// Post receives first, then pack/post the sends.  The group buffers belong to
// the plan and stay valid until finishRecvs()/finishPlan().
template <typename Plan, typename PackFn>
void postPlan(Plan& plan, std::vector<MPI_Request>& sendRequests,
              std::vector<MPI_Request>& recvRequests, PackFn&& pack) {
  postRecvs(plan, recvRequests);
  postSends(plan, sendRequests, pack);
}

// Wait for every send of a posted plan.  Kept separate so the hidden step can
// wait after the interior collision and before the shared-memory halo copy.
inline void waitSends(std::vector<MPI_Request>& sendRequests) {
  if (!sendRequests.empty()) {
    MPI_Waitall(static_cast<int>(sendRequests.size()), sendRequests.data(),
                MPI_STATUSES_IGNORE);
  }
}

// Wait for the receives and unpack them.  For small plans, keep the old
// per-group wait so the first group can be unpacked before the last arrives;
// for larger plans, wait all and unpack the groups in parallel.  MPI waits
// stay outside the parallel region.
template <typename Plan, typename UnpackFn>
void finishRecvs(Plan& plan, std::vector<MPI_Request>& recvRequests, UnpackFn&& unpack) {
  const int groups = static_cast<int>(plan.Recv.size());
  std::size_t total = 0;
  for (const auto& group : plan.Recv) total += group.Buffer.size();
  const bool parallel = groups >= kParallelMinGroups && total >= kParallelMinElements;
  if (parallel) {
    if (!recvRequests.empty()) {
      MPI_Waitall(static_cast<int>(recvRequests.size()), recvRequests.data(),
                  MPI_STATUSES_IGNORE);
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
    for (int g = 0; g < groups; ++g) {
      const auto& group = plan.Recv[g];
      for (const CommRef& ref : group.Refs) {
        unpack(ref, group.Buffer, ref.Offset);
      }
    }
  } else {
    for (int g = 0; g < groups; ++g) {
      MPI_Wait(&recvRequests[g], MPI_STATUS_IGNORE);
      const auto& group = plan.Recv[g];
      for (const CommRef& ref : group.Refs) {
        unpack(ref, group.Buffer, ref.Offset);
      }
    }
  }
  recvRequests.clear();
}

// The synchronous form: post, wait for the sends, wait/unpack the receives.
template <typename Plan, typename PackFn, typename UnpackFn>
void runPlan(Plan& plan, PackFn&& pack, UnpackFn&& unpack) {
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> recvRequests;
  postPlan(plan, sendRequests, recvRequests, pack);
  waitSends(sendRequests);
  finishRecvs(plan, recvRequests, unpack);
}

}  // namespace flbcomm

#endif  // MPI_ENABLED
