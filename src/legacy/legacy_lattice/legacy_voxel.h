
#pragma once

#include "data_struct/Vector.h"

// single voxel(3D), can be used as pixel(2D)
template <typename T, unsigned int D>
class Voxel : public Vector<T, D> {
 private:
  // index in std::vector<Voxel<T, 3>> _Voxels
  int id_;
  // voexl flag
  std::uint8_t flag_;
  // neighbor
  std::vector<Voxel<T, D> *> neighbor_;
  // neighbor count
  int neighbor_count_;

 public:
  // constructors
  Voxel() : Vector<T, D>(), id_(-1), flag_(0), neighbor_count_(0) {}
  // Voxel(T value, int id = 0, int flag = 0)
  //     : Vector<T, D>(value), id_(id), flag_(flag), neighbor_count_(0) {}
  Voxel(const Vector<T, D> &vec, int id = 0, std::uint8_t flag = 0)
      : Vector<T, D>(vec), id_(id), flag_(flag), neighbor_count_(0) {}
  template <typename... Args>
  Voxel(int id, std::uint8_t flag, Args... args)
      : Vector<T, D>(args...), id_(id), flag_(flag), neighbor_count_(0) {}

  // get index in std::vector<Voxel<T, d>> _Voxels
  const int &getId() const { return id_; }
  inline const Vector<T, D> &getLoc() const {
    return static_cast<const Vector<T, D> &>(*this);
  }
  inline const int &getFlag() const { return flag_; }
  // set
  void setFlag(std::uint8_t flag) { flag_ = flag; }
  // get neighbor container
  const std::vector<Voxel<T, D> *> &getNeighborList() const {
    return neighbor_;
  }
  // get neighbor voxel
  Voxel<T, D> *getNeighbor(int i) { return neighbor_[i]; }
  const Voxel<T, D> *getNeighbor(int i) const { return neighbor_[i]; }
  int getNeighborId(int i) const { return neighbor_[i]->getId(); }
  // get neighbor count
  int getNeighborCount() const { return neighbor_count_; }
  // set neighbor voxel
  void initNeighborList(int size) { neighbor_.resize(size, nullptr); }
  // add is not use for now cause initNeighborList is more convenient
  // void addNeighbor(Voxel<T, D> *voxel) {
  //   neighbor_.push_back(voxel);
  //   neighbor_count_++;
  // }
  void setNeighbor(int i, Voxel<T, D> *voxel) {
    if (neighbor_[i] == nullptr) neighbor_count_++;
    neighbor_[i] = voxel;
  }

  // get index
  // get index in a cube with size (1, Nx, Nx * Ny) or (1, Nx)
  // int getIdx(const Vector<int, D> &size) const {
  //   int idx = 0;
  //   for (unsigned int i = 0; i < D; ++i) {
  //     idx += this->operator[](i) * size[i];
  //   }
  //   return idx;
  // }
};