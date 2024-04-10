// legacy_lattice_set.h

#pragma once

#include "lbm/lattice_set.h"

namespace lat {

template <typename LatSet>
struct LatStru_Set2D {
  static inline void Set_Nbr(int ni, int *Nbr) {
    for (int i = 0; i < LatSet::q; i++) {
      Nbr[i] = Index2D::GetId(LatSet::c[i][0], LatSet::c[i][1], ni);
    }
  }
  static void Set_Period(int ni, int *Per);
};

template <typename LatSet>
struct LatStru_Set3D {
  static inline void Set_Nbr(int ni, int ninj, int *Nbr) {
    for (int i = 0; i < LatSet::q; i++) {
      Nbr[i] = Index3D::GetId(LatSet::c[i][0], LatSet::c[i][1], LatSet::c[i][2],
                              ni, ninj);
    }
  }
};

template <typename LatSet>
struct D1Q3 : public Basic_Lattice_Set<1, 3> {
  int Ni;
  int Nbr[3];   // neighbor
  int Per[3];   // periodic
  int VPer[3];  // generalised periodic, use virtual nodes
  int Sym[3];   // symmetry
  // constructor
  D1Q3(int ni) : Ni(ni), Nbr{0, 1, -1} {}
  void SetPeriod() {
    Per[0] = 0;
    Per[1] = Ni - 3;
    Per[2] = 3 - Ni;
  }
  void SetGeneralPeriod() {
    VPer[0] = 0;
    VPer[1] = Ni - 4;
    VPer[2] = 4 - Ni;
  }
  void SetSymmetry() {
    // transversal axis of symmetry
    Sym[1] = 2;
    Sym[2] = 1;
  }
};

template <typename LatSet>
struct D2Q4 : public Basic_Lattice_Set<2, 4> {
  int Ni;
  int Nj;
  int Nbr[q];

  // constructor
  // D2Q4(int ni, int nj) : Ni(ni), Nj(nj), Nbr{1, Ni, -1, -Ni} {}
  D2Q4(int ni, int nj) : Ni(ni), Nj(nj) {
    LatStru_Set2D<LatSet>::Set_Nbr(ni, Nbr);
  }
};

struct D2Q5 : public Basic_Lattice_Set<2, 5> {
  int Ni;
  int Nj;
  int Nbr[5];   // neighbor
  int Per[5];   // periodic
  int VPer[5];  // generalised periodic, use virtual nodes
  int Sym[5];   // symmetry
  // constructor
  D2Q5(int ni, int nj) : Ni(ni), Nj(nj), Nbr{0, 1, Ni, -1, -Ni} {}
  void SetPeriod() {
    Per[0] = 0;
    Per[1] = Ni - 3;
    Per[2] = Ni * (Nj - 3);
    Per[3] = 3 - Ni;
    Per[4] = -Ni * (Nj - 3);
  }
  void SetGeneralPeriod() {
    VPer[0] = 0;
    VPer[1] = Ni - 4;
    VPer[2] = Ni * (Nj - 4);
    VPer[3] = 4 - Ni;
    VPer[4] = -Ni * (Nj - 4);
  }
  void SetSymmetry() {
    // transversal axis of symmetry 2-4
    Sym[2] = 4;
    Sym[4] = 2;
    // longitudinal axis of symmetry 1-3
    Sym[1] = 3;
    Sym[3] = 1;
  }
};

struct D2Q9 : public Basic_Lattice_Set<2, 9> {
  int Ni;
  int Nj;
  int Nbr[9];   // neighbor
  int Per[9];   // periodic
  int VPer[9];  // generalised periodic, use virtual nodes
  int Sym[9];   // symmetry
  // constructor
  D2Q9(int ni, int nj)
      : Ni(ni),
        Nj(nj),
        Nbr{0, 1, Ni, -1, -Ni, Ni + 1, Ni - 1, -Ni - 1, -Ni + 1} {}
  // 0: left and right, 1: top and bottom
  void SetPeriod(int type) {
    Per[0] = 0;
    Per[1] = Ni - 3;
    Per[2] = Ni * (Nj - 3);
    Per[3] = 3 - Ni;
    Per[4] = -Ni * (Nj - 3);
    // left and right
    if (type == 0) {
      Per[5] = Per[1] - Ni;
      Per[8] = Per[1] + Ni;
      Per[6] = Per[3] - Ni;
      Per[7] = Per[3] + Ni;
    }
    // top and bottom
    if (type == 1) {
      Per[7] = Per[4] + 1;
      Per[8] = Per[4] - 1;
      Per[5] = Per[2] - 1;
      Per[6] = Per[2] + 1;
    }
  }
  // 0: left and right, 1: top and bottom
  void SetGeneralPeriod(int type) {
    VPer[0] = 0;
    VPer[1] = Ni - 4;
    VPer[2] = Ni * (Nj - 4);
    VPer[3] = 4 - Ni;
    VPer[4] = -Ni * (Nj - 4);
    // left and right
    if (type == 0) {
      VPer[5] = VPer[1] - Ni;
      VPer[8] = VPer[1] + Ni;
      VPer[6] = VPer[3] - Ni;
      VPer[7] = VPer[3] + Ni;
    }
    // top and bottom
    if (type == 1) {
      VPer[7] = VPer[4] + 1;
      VPer[8] = VPer[4] - 1;
      VPer[5] = VPer[2] - 1;
      VPer[6] = VPer[2] + 1;
    }
  }
  void SetSymmetry(int type) {
    // transversal axis of symmetry
    if (type == 0) {
      // 2-4 5-8 6-7
      Sym[2] = 4;
      Sym[4] = 2;
      Sym[5] = 8;
      Sym[6] = 7;
      Sym[7] = 6;
      Sym[8] = 5;
    }
    // longitudinal axis of symmetry
    if (type == 1) {
      // 1-3 5-6 7-8
      Sym[1] = 3;
      Sym[3] = 1;
      Sym[5] = 6;
      Sym[6] = 5;
      Sym[7] = 8;
      Sym[8] = 7;
    }
  }
};

template <typename LatSet>
struct D3Q7 : public Basic_Lattice_Set<3, 7> {
  int Ni;
  int Nj;
  int Nk;
  int Nbr[7];   // neighbor
  int Per[7];   // periodic
  int VPer[7];  // generalised periodic, use virtual nodes
  int Sym[7];   // symmetry
  D3Q7(int ni, int nj, int nk) : Ni(ni), Nj(nj), Nk(nk) {
    LatStru_Set3D<LatSet>::Set_Nbr(ni, ni * nj, Nbr);
  }
};

template <typename LatSet>
struct D3Q15 : public Basic_Lattice_Set<3, 15> {
  int Ni;
  int Nj;
  int Nk;
  int Nbr[15];  // neighbor

  // constructor
  D3Q15(int ni, int nj, int nk) : Ni(ni), Nj(nj), Nk(nk) {
    LatStru_Set3D<LatSet>::Set_Nbr(ni, ni * nj, Nbr);
  }
};

template <typename LatSet>
struct D3Q27 : public Basic_Lattice_Set<3, 27> {
  int Ni;
  int Nj;
  int Nk;
  int Nbr[27];  // neighbor

  // constructor
  D3Q27(int ni, int nj, int nk) : Ni(ni), Nj(nj), Nk(nk) {
    LatStru_Set3D<LatSet>::Set_Nbr(ni, ni * nj, Nbr);
  }
};

}  // namespace lat

namespace ca {

template <int D, int Q>
struct Basic_Lattice_Set {
  // DdQq
  static constexpr int d = D;
  static constexpr int q = Q;
};

// Von Neumann neighborhood type
// D2Q4
struct D2Q4 : public Basic_Lattice_Set<2, 4> {
  static constexpr int c[q][d] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
  static constexpr int opp[q] = {2, 3, 0, 1};
};

// Moore neighborhood type
// D2Q8
struct D2Q8 : public Basic_Lattice_Set<2, 8> {
  static constexpr int c[q][d] = {{1, 0}, {0, 1},  {-1, 0},  {0, -1},
                                  {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  static constexpr int opp[q] = {2, 3, 0, 1, 6, 7, 4, 5};
};

}  // namespace ca

namespace calat {

struct D2Q4 : public Basic_Lattice_Set<2, 4> {
  int Ni;
  int Nj;
  int Nbr[4];  // neighbor

  // constructor
  D2Q4(int ni, int nj) : Ni(ni), Nj(nj), Nbr{1, Ni, -1, -Ni} {}
};

struct D2Q8 : public Basic_Lattice_Set<2, 8> {
  int Ni;
  int Nj;
  int Nbr[8];  // neighbor

  // constructor
  D2Q8(int ni, int nj)
      : Ni(ni), Nj(nj), Nbr{1, Ni, -1, -Ni, 1 + Ni, -1 + Ni, -1 - Ni, 1 - Ni} {}
};

}  // namespace calat