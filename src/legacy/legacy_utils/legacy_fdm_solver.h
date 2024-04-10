#pragma once
#include "legacy/legacy_lattice/legacy_voxel.h"
// #include "legacy/lattice/legacyfield3D.h"
#include "util.h"

namespace FDM {
template <typename T>
class FDM3D {
 private:
  std::vector<Voxel<T, 3>> &_Voxels;
  T *f;

 public:
  FDM3D(std::vector<Voxel<T, 3>> &voxels, T *f_) : _Voxels(voxels), f(f_) {}
  // FDM parser
  // first order derivative: c[q] = {0, 0, 0},  {1, 0, 0}, {-1, 0, 0}, {0, 1,
  // 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}  0-6
  // partial x, central difference
  inline T p_x(int id) const {
    const Voxel<T, 3> &voxel = _Voxels[id];
    // check if voxel is boundary cell
    if (voxel.getFlag() == 0) {
      // if not, use central difference
      return (f[voxel.getNeighbor(1)->getId()] -
              f[voxel.getNeighbor(2)->getId()]) *
             T(0.5);
    } else {
      // if yes, use forward/ backward difference
      if (voxel.getNeighbor(1) == nullptr && voxel.getNeighbor(2) != nullptr)
        return f[id] - f[voxel.getNeighbor(2)->getId()];
      else if (voxel.getNeighbor(1) != nullptr &&
               voxel.getNeighbor(2) == nullptr)
        return f[voxel.getNeighbor(1)->getId()] - f[id];
      else if (voxel.getNeighbor(1) != nullptr &&
               voxel.getNeighbor(2) != nullptr)
        return (f[voxel.getNeighbor(1)->getId()] -
                f[voxel.getNeighbor(2)->getId()]) *
               T(0.5);
      else
        return 0;
    }
  }
  // partial y, central difference
  inline T p_y(int id) const {
    const Voxel<T, 3> &voxel = _Voxels[id];
    // check if voxel is boundary cell
    if (voxel.getFlag() == 0) {
      // if not, use central difference
      return (f[voxel.getNeighbor(3)->getId()] -
              f[voxel.getNeighbor(4)->getId()]) *
             T(0.5);
    } else {
      // if yes, use forward/ backward difference
      if (voxel.getNeighbor(3) == nullptr && voxel.getNeighbor(4) != nullptr)
        return f[id] - f[voxel.getNeighbor(4)->getId()];
      else if (voxel.getNeighbor(3) != nullptr &&
               voxel.getNeighbor(4) == nullptr)
        return f[voxel.getNeighbor(3)->getId()] - f[id];
      else if (voxel.getNeighbor(3) != nullptr &&
               voxel.getNeighbor(4) != nullptr)
        return (f[voxel.getNeighbor(3)->getId()] -
                f[voxel.getNeighbor(4)->getId()]) *
               T(0.5);
      else
        return 0;
    }
  }
  // partial z, central difference
  inline T p_z(int id) const {
    const Voxel<T, 3> &voxel = _Voxels[id];
    // check if voxel is boundary cell
    if (voxel.getFlag() == 0) {
      // if not, use central difference
      return (f[voxel.getNeighbor(5)->getId()] -
              f[voxel.getNeighbor(6)->getId()]) *
             T(0.5);
    } else {
      // if yes, use forward/ backward difference
      if (voxel.getNeighbor(5) == nullptr && voxel.getNeighbor(6) != nullptr)
        return f[id] - f[voxel.getNeighbor(6)->getId()];
      else if (voxel.getNeighbor(5) != nullptr &&
               voxel.getNeighbor(6) == nullptr)
        return f[voxel.getNeighbor(5)->getId()] - f[id];
      else if (voxel.getNeighbor(5) != nullptr &&
               voxel.getNeighbor(6) != nullptr)
        return (f[voxel.getNeighbor(5)->getId()] -
                f[voxel.getNeighbor(6)->getId()]) *
               T(0.5);
      else
        return 0;
    }
  }
  // get gradient: grad = {partial x, partial y, partial z}
  inline Vector<T, 3> grad(int id) const {
    return Vector<T, 3>(p_x(id), p_y(id), p_z(id));
  }
  // get Normalizedgradient = grad / |grad|
  inline Vector<T, 3> ngrad(int id) const {
    Vector<T, 3> grad = grad(id);
    return grad / grad.norm();
  }
  // partial xx
  inline T p_xx(int id) const {
    const Voxel<T, 3> &voxel = _Voxels[id];
    // check if voxel is boundary cell
    if (voxel.getFlag() == 0) {
      return (f[voxel.getNeighbor(1)->getId()] +
              f[voxel.getNeighbor(2)->getId()] - 2 * f[id]);
    } else {
      if (voxel.getNeighbor(1) != nullptr && voxel.getNeighbor(2) != nullptr)
        return (f[voxel.getNeighbor(1)->getId()] +
                f[voxel.getNeighbor(2)->getId()] - 2 * f[id]);
      else
        return 0;
    }
  }
  // partial yy
  inline T p_yy(int id) const {
    const Voxel<T, 3> &voxel = _Voxels[id];
    // check if voxel is boundary cell
    if (voxel.getFlag() == 0) {
      return (f[voxel.getNeighbor(3)->getId()] +
              f[voxel.getNeighbor(4)->getId()] - 2 * f[id]);
    } else {
      if (voxel.getNeighbor(3) != nullptr && voxel.getNeighbor(4) != nullptr)
        return (f[voxel.getNeighbor(3)->getId()] +
                f[voxel.getNeighbor(4)->getId()] - 2 * f[id]);
      else
        return 0;
    }
  }
  // partial zz
  inline T p_zz(int id) const {
    const Voxel<T, 3> &voxel = _Voxels[id];
    // check if voxel is boundary cell
    if (voxel.getFlag() == 0) {
      return (f[voxel.getNeighbor(5)->getId()] +
              f[voxel.getNeighbor(6)->getId()] - 2 * f[id]);
    } else {
      if (voxel.getNeighbor(5) != nullptr && voxel.getNeighbor(5) != nullptr)
        return (f[voxel.getNeighbor(5)->getId()] +
                f[voxel.getNeighbor(6)->getId()] - 2 * f[id]);
      else
        return 0;
    }
  }
  // get Inabla of gradient
  // Inabla = (partial_x, 0, 0; 0, partial_y, 0; 0, 0, partial_z)
  // Inabla(grad) = (partial_xx, partial_xy, partial_xz)^T
  Vector<T, 3> Inablagrad(int id) const {
    return Vector<T, 3>(p_xx(id), p_yy(id), p_zz(id));
  }
  // get divergence of gradient
  // divgrad = partial_xx + partial_yy + partial_zz
  T divgrad(int id) const { return p_xx(id) + p_yy(id) + p_zz(id); }
  // {1, 1, 0},  {-1, -1, 0},  {1, 0, 1},  {-1, 0, -1},
  // {0, 1, 1},  {0, -1, -1},  {1, -1, 0}, {-1, 1, 0},
  // {1, 0, -1}, {-1, 0, 1},   {0, 1, -1}, {0, -1, 1}, 7-18
  // partial xy
  inline T p_xy(int id) const {
    const Voxel<T, 3> &voxel = _Voxels[id];
    // check if voxel is boundary cell
    if (voxel.getFlag() == 0) {
      return (f[voxel.getNeighbor(7)->getId()] +
              f[voxel.getNeighbor(8)->getId()] -
              f[voxel.getNeighbor(13)->getId()] -
              f[voxel.getNeighbor(14)->getId()]) *
             T(0.25);
    } else {
      if (voxel.getNeighbor(7) != nullptr && voxel.getNeighbor(8) != nullptr &&
          voxel.getNeighbor(13) != nullptr && voxel.getNeighbor(14) != nullptr)
        return (f[voxel.getNeighbor(7)->getId()] +
                f[voxel.getNeighbor(8)->getId()] -
                f[voxel.getNeighbor(13)->getId()] -
                f[voxel.getNeighbor(14)->getId()]) *
               T(0.25);
      else
        return 0;
    }
  }
  inline T p_xz(int id) const {
    const Voxel<T, 3> &voxel = _Voxels[id];
    // check if voxel is boundary cell
    if (voxel.getFlag() == 0) {
      return (f[voxel.getNeighbor(9)->getId()] +
              f[voxel.getNeighbor(10)->getId()] -
              f[voxel.getNeighbor(15)->getId()] -
              f[voxel.getNeighbor(16)->getId()]) *
             T(0.25);
    } else {
      if (voxel.getNeighbor(9) != nullptr && voxel.getNeighbor(10) != nullptr &&
          voxel.getNeighbor(15) != nullptr && voxel.getNeighbor(16) != nullptr)
        return (f[voxel.getNeighbor(9)->getId()] +
                f[voxel.getNeighbor(10)->getId()] -
                f[voxel.getNeighbor(15)->getId()] -
                f[voxel.getNeighbor(16)->getId()]) *
               T(0.25);
      else
        return 0;
    }
  }
  inline T p_yz(int id) const {
    const Voxel<T, 3> &voxel = _Voxels[id];
    // check if voxel is boundary cell
    if (voxel.getFlag() == 0) {
      return (f[voxel.getNeighbor(11)->getId()] +
              f[voxel.getNeighbor(12)->getId()] -
              f[voxel.getNeighbor(17)->getId()] -
              f[voxel.getNeighbor(18)->getId()]) *
             T(0.25);
    } else {
      if (voxel.getNeighbor(11) != nullptr &&
          voxel.getNeighbor(12) != nullptr &&
          voxel.getNeighbor(17) != nullptr && voxel.getNeighbor(18) != nullptr)
        return (f[voxel.getNeighbor(11)->getId()] +
                f[voxel.getNeighbor(12)->getId()] -
                f[voxel.getNeighbor(17)->getId()] -
                f[voxel.getNeighbor(18)->getId()]) *
               T(0.25);
      else
        return 0;
    }
  }
};

}  // namespace FDM
