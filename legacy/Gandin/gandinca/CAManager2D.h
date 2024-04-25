// CAManager2D.h
#pragma once

#include "grow.h"
#include "legacy/legacy_lbm/legacy_lbmanager2d.h"
#include "nuc.h"

template <typename T, typename LatStru>
class CAManager2D : public Nuc2D<T, LatStru>, public Grow2D<T, LatStru> {
 private:
  GetGCells2D& Rem;
  LBManager2D<T>& LBM;

 public:
  CAManager2D(CAGField2D<T>& ca, GandinConverter<T>& convca, LBManager2D<T>& lbm,
              GetGCells2D& remains)
      : Rem(remains),
        LBM(lbm),
        Nuc2D<T, LatStru>(convca, lbm, ca),
        Grow2D<T, LatStru>(ca, convca, lbm) {}

  // remember to call Get_remain cells before lbm
  void apply() {
    
    // Rem.Get();
    Nuc2D<T, LatStru>::Create_PreNucs(Rem.Get_Bulks(), Rem.Get_Surfs());
    Nuc2D<T, LatStru>::Nucleation(Grow2D<T, LatStru>::getNewNucs(), LBM);
    Grow2D<T, LatStru>::Grow();
    Grow2D<T, LatStru>::Capture();
    Grow2D<T, LatStru>::addto_NewGrowings();
  }

  // remember to call Get_remain cells  before lbm
  void SingleNuc(int latx, int laty,
                 T orine = (std::rand() % 90) * M_PI / T(180)) {
    Nuc2D<T, LatStru>::Single_Nucleation(latx, laty, orine,
                                         Grow2D<T, LatStru>::getNewNucs());
  }

  void apply_SingleNuc() {
    // Rem.Get();
    Grow2D<T, LatStru>::Grow();
    Grow2D<T, LatStru>::Capture();
    Grow2D<T, LatStru>::addto_NewGrowings();
  }
};