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

// block_link_meta.h

#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "geometry/basic_geometry.h"
#include "parallel/communicator.h"
#include "parallel/link_key.h"

namespace flbgeo {

// The three kinds of same-level link.  The receiver's block decides the split:
// a face comm excludes corner/edge cells, a corner comm is one cell, and a 3D
// edge comm is the edge group.  The sender uses the receiver's partition, which
// is why the same builder is called with the receiver as `self` on both sides.
enum class LinkKind : std::uint8_t { Face = 0, Corner = 1, Edge = 2 };

// One directional link between two same-level blocks, described from the
// receiver's side.  `Recv` indexes the receiver's ghost layer, `PeerSend`
// indexes the sender's base layer; the two lists are aligned pair by pair.
// `Order` is derived only from the receiver/sender block ids and the receiver's
// link kind/index, so both ends derive the same value without exchanging it.
struct SameLevelLink {
  LinkKind Kind{LinkKind::Face};
  std::uint16_t Index{};
  NbrDirection Direction{NbrDirection::NONE};
  flbcomm::LinkKey Order{};
  std::vector<std::size_t> Recv;
  std::vector<std::size_t> PeerSend;
};

// Both ends call this with the receiver as `receiver`, so the key is identical
// without exchanging it.  The four fields keep the order exact: two distinct
// links can never compare equal, unlike a hash.
inline flbcomm::LinkKey linkOrder(int receiver, int sender, LinkKind kind,
                                  std::uint16_t index) {
  return flbcomm::LinkKey{receiver, sender, static_cast<std::uint16_t>(kind),
                          index};
}

// Enumerate the directional links from `self` (receiver) to `nbr` (sender).
// The receiver's cell lists are produced by self/self's extended block; the
// sender's send list is produced in the sender's block coordinates.  Both ends
// call this with the receiver as `self`, so the partition, direction and Order
// are identical on the two ranks.
template <typename T, unsigned int D>
std::vector<SameLevelLink> BuildSameLevelLinks(const BasicBlock<T, D>& selfBase,
                                               int selfOverlap,
                                               const BasicBlock<T, D>& nbrBase,
                                               int nbrOverlap) {
  std::vector<SameLevelLink> links;

  const BasicBlock<T, D> selfExt = selfBase.getExtBlock(selfOverlap);
  const BasicBlock<T, D> nbrExt = nbrBase.getExtBlock(nbrOverlap);
  const BasicBlock<T, D> selfExt1 = selfBase.getExtBlock(1);

  std::vector<std::size_t> recv;
  std::vector<std::size_t> peerSend;
  selfExt.getCellIdx(selfExt1, nbrBase, recv);
  nbrExt.getCellIdx(nbrBase, selfExt1, peerSend);

  std::vector<std::size_t> cornerRecv;
  std::vector<std::size_t> cornerSend;
  selfExt.ExcludeCornerIdx(recv, peerSend, cornerRecv, cornerSend);

  std::vector<std::vector<std::size_t>> edgeRecv;
  std::vector<std::vector<std::size_t>> edgeSend;
  if constexpr (D == 3) {
    selfExt.ExcludeEdgeIdx(recv, peerSend, edgeRecv, edgeSend);
  }

  const int receiver = selfBase.getBlockId();
  const int sender = nbrBase.getBlockId();

  if (!recv.empty()) {
    SameLevelLink link;
    link.Kind = LinkKind::Face;
    link.Index = 0;
    link.Order = linkOrder(receiver, sender, link.Kind, link.Index);
    if constexpr (D == 2) {
      link.Direction = getEdgeNbrDirection<2>(selfExt.whichEdge(recv[0]));
    } else {
      link.Direction = getFaceNbrDirection(selfExt.whichFace(recv[0]));
    }
    link.Recv = std::move(recv);
    link.PeerSend = std::move(peerSend);
    links.push_back(std::move(link));
  }

  for (std::size_t i = 0; i < cornerRecv.size(); ++i) {
    const int corner = selfExt.whichCorner(cornerRecv[i]);
    SameLevelLink link;
    link.Kind = LinkKind::Corner;
    link.Index = static_cast<std::uint16_t>(corner);
    link.Order = linkOrder(receiver, sender, link.Kind, link.Index);
    link.Direction = getCornerNbrDirection<D>(corner);
    link.Recv.push_back(cornerRecv[i]);
    link.PeerSend.push_back(cornerSend[i]);
    links.push_back(std::move(link));
  }

  if constexpr (D == 3) {
    for (int edge = 0; edge < static_cast<int>(edgeRecv.size()); ++edge) {
      if (edgeRecv[edge].empty()) continue;
      SameLevelLink link;
      link.Kind = LinkKind::Edge;
      link.Index = static_cast<std::uint16_t>(edge);
      link.Order = linkOrder(receiver, sender, link.Kind, link.Index);
      link.Direction = getEdgeNbrDirection<3>(selfExt.whichEdge(edgeRecv[edge][0]));
      link.Recv = std::move(edgeRecv[edge]);
      link.PeerSend = std::move(edgeSend[edge]);
      links.push_back(std::move(link));
    }
  }

  return links;
}

}  // namespace flbgeo
