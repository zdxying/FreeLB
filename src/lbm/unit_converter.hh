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

#pragma once

#include "lbm/unit_converter.h"

template <typename T>
void BaseConverter<T>::check(int &check_status) {
  std::cout << "[BasicConverter]:\n"
            << "deltaX:           " << deltaX << "\n"
            << "deltaT:           " << deltaT << "\n"
            << "charL:            " << charL << "\n"
            << "charU:            " << charU << "\n"
            << "VisKine:          " << VisKine << "\n"
            << "Re:               " << Re << "\n"
            << "Conv_Acc:         " << Conv_Acc << "\n"
            << "Lattice_charU:    " << Lattice_charU << "\n"
            << "Lattice_charL:    " << Lattice_charL << "\n"
            << "Lattice_VisKine:  " << Lattice_VisKine << "\n"
            << "Lattice_Re:       " << Lattice_Re << "\n";
  if (Lattice_Re > 100) {
    std::cout << "[Warning]: Lattice_Re > 100, Press any key to continue..."
              << std::endl;
    getchar();
    check_status = 1;
  }

  if (Lattice_charU > 0.4) {
    std::cout << "Error: Lattice_charU > 0.4" << std::endl;
    exit(-1);
  } else {
    if (Lattice_RT < 0.55) {
      T Minstable_Lattice_RT =
          0.5 + 0.125 * Lattice_charU;  // 0.5 + 0.125 * 0.4 = 0.55, 0.5 + 0.125
                                        // * 0.03 = 0.50375
      if (Lattice_RT < Minstable_Lattice_RT) {
        std::cout << "Error: Lattice_charU too large, Minstable_Lattice_RT = "
                  << Minstable_Lattice_RT << std::endl;
        exit(-1);
      }
    }
  }
}

template <typename T>
void TempConverter<T>::check(int &check_status) {
  std::cout << "[TempConverter]:\n"
            << "Tl:               " << Tl << "\n"
            << "Th:               " << Th << "\n"
            << "Conv_dT:          " << Conv_dT << "\n"
            << "TInit:            " << TInit << "\n"
            << "Lattice_TInit:    " << Lattice_TInit << "\n";
  if (Lattice_RTT < 0.5) {
    std::cout << "Error: Lattice_RTT < 0.5" << std::endl;
    exit(-1);
  }

  if (SHeatCapConverted) {
    Lattice_SHeatCap = static_cast<T>(1e9) * SHeatCap / BaseConv.Conv_U /
                       BaseConv.Conv_U * Conv_dT;
    std::cout << "Lattice_SHeatCap: " << Lattice_SHeatCap << "\n";
  }
  if (LatentHeatConverted) {
    Lattice_LatentHeat =
        static_cast<T>(1e9) * LatentHeat / BaseConv.Conv_U / BaseConv.Conv_U;
    std::cout << "Lattice_LatHeat:  " << Lattice_LatentHeat << "\n";
  }
  if (SHeatCapConverted && LatentHeatConverted) {
    Lattice_LatHeat_SHeatCap = Lattice_LatentHeat / Lattice_SHeatCap;
    std::cout << "Lat_LatHeat_SHC:  " << Lattice_LatHeat_SHeatCap << "\n";
  }
  if (TExpanConverted) {
    std::cout << "Lattice_gbetaT:   " << Lattice_gbetaT << "\n"
              << "Ra:               " << Ra << "\n";
    if (Lattice_gbetaT > T(0.001) || Lattice_gbetaT < T(-0.001)) {
      std::cout << "Lattice_gbetaT too large, Press any key to continue... "
                << std::endl;
      getchar();
      check_status = 1;
    }
  }
}

template <typename T>
void ConcConverter<T>::check(int &check_status) {
  std::cout << "[ConcConverter]:\n"
            << "Cl:               " << Cl << "\n"
            << "Ch:               " << Ch << "\n"
            << "Conv_dC:          " << Conv_dC << "\n"
            << "CInit:            " << CInit << "\n"
            << "Lattice_CInit:    " << Lattice_CInit << "\n";

  if (Lattice_RTC < T(0.5)) {
    std::cout << "Error: Lattice_RTC = " << Lattice_RTC << " < 0.5"
              << std::endl;
    exit(-1);
  }
  if (CExpanConverted) {
    std::cout << "Lattice_gbetaC:   " << Lattice_gbetaC << "\n";

    if (Lattice_gbetaC > T(0.01) || Lattice_gbetaC < T(-0.01)) {
      std::cout << "Lattice_gbetaC too large, Press any key to continue..."
                << std::endl;
      getchar();
      check_status = 1;
    }
  }
}

template <typename T>
void ZSConverter<T>::check(int &check_status) {
  std::cout << "[ZS Model Converter]:\n"
            << "Lattice_GT_Coef: " << Lattice_GT_Coef << "\n";
}

template <typename T>
void GandinConverter<T>::check(int &check_status) {
  std::cout << "[GandinConverter]:  \n"
            << "Lat_DT_Mean_Bulk: " << Lattice_DT_Mean_Bulk << "\n"
            << "Lat_DT_Std_Bulk:  " << Lattice_DT_Std_Bulk << "\n"
            << "Lat_DT_Mean_Surf: " << Lattice_DT_Mean_Surf << "\n"
            << "Lat_DT_Std_Surf:  " << Lattice_DT_Std_Surf << "\n"
            << "Lat_Nuc_Dens_Bulk:" << Lattice_NucDens_Bulk << "\n"
            << "Lat_Nuc_Dens_Surf:" << Lattice_NucDens_Surf << "\n"
            << "Lat_GrowthPara:   " << Lattice_GrowthPara << "\n";
  T MaxGrowth = Lattice_GrowthPara *
                (Lattice_DT_Mean_Bulk + 3 * Lattice_DT_Std_Bulk) *
                (Lattice_DT_Mean_Bulk + 3 * Lattice_DT_Std_Bulk);
  if (Lattice_NucDens_Bulk >= 1) {
    std::cout << "Error: Lattice_NucDens_Bulk >= 1" << std::endl;
    exit(-1);
  }
  if (Lattice_NucDens_Surf >= 1) {
    std::cout << "Error: Lattice_NucDens_Surf >= 1" << std::endl;
    exit(-1);
  }
  if (MaxGrowth > 0.2) {
    std::cout << "Error: MaxGrowth too large = " << MaxGrowth << std::endl;
    exit(-1);
  }
  if (TempConv.SHeatCapConverted && TempConv.LatentHeatConverted) {
    T MaxSourceT = MaxGrowth * TempConv.Lattice_LatHeat_SHeatCap;
    if (MaxSourceT > 0.1) {
      std::cout << "Error: MaxSourceT too large = " << MaxSourceT << std::endl;
      exit(-1);
    }
  }
}

template <typename T>
void UnitConvManager<T>::Check_and_Print() {
  MPI_RANK(0)
  int check_status = 0;
  Printer::Print_Banner("Convert Log");

  BaseConv->check(check_status);
  if (TempConv != nullptr) {
    TempConv->check(check_status);
  }
  if (ConcConv != nullptr) {
    ConcConv->check(check_status);
  }
  if (CAConv != nullptr) {
    CAConv->check(check_status);
  }

  if (check_status == 0) {
    std::cout << "Simulation parameters correctly set." << std::endl;
  } else {
    std::cout << "Simulation parameters set with warnings." << std::endl;
  }

  std::cout << "[Lattice Parameters]:\n"
            << "Lattice_RT =     " << BaseConv->getLattice_RT() << "\n"
            << "OMEGA =          " << BaseConv->getOMEGA() << "\n";
  if (TempConv != nullptr) {
    std::cout << "Lattice_RTT =    " << TempConv->getLattice_RT() << "\n"
              << "OMEGAT =         " << TempConv->getOMEGA() << "\n";
  }
  if (ConcConv != nullptr) {
    std::cout << "Lattice_RTC =    " << ConcConv->getLattice_RT() << "\n"
              << "OMEGAC =         " << ConcConv->getOMEGA() << "\n";
  }
  std::cout << "-------------------------------------------------\n"
            << std::endl;
}