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

/*Lattice Boltzmann Unit Converter*/
#pragma once

#include <iostream>

#include "data_struct/Vector.h"

// Lattice_deltaX = 1   ->  Conv_L = deltaX
// Lattice_deltaT = 1   ->  Conv_Time = deltaT
// Lattice_rho = 1   ->  Conv_rho = rho
template <typename T>
struct AbstractConverter {
  // /*----------------------------phys param----------------------------*/
  // /*-----------*/
  // T deltaX;  // mm          // spacing between two lattice cells
  // T deltaT;  // s           // time step
  // T charL;   // mm          // i.e. simulation domain size
  // T charU;   // mm / s      // maximal or expected velocity during simulation
  // T charP;   // g / (mm s^2) = Pa
  // T VisKine; // mm^2 / s  // kinematic viscosity
  // T rho;     // g / mm^3    // density
  // T Re;      // Reynolds number, Re = U0*L0/nu = charU*charL/VisKine
  // /*-----------*/
  // T Th;               // K      //characteristic physical high temperature
  // T Tl;               // K      //characteristic physical low temperature
  // T TCond;            // W / (mm K)// thermal conductivity
  // T SHeatCap;         // J / (g K) = 10^9 * mm^2 / s^2 / K // specific heat
  // capacity at constant pressure T TDiff;            // mm^2 / s  // thermal
  // diffusion T Ra;               // Rayleigh number, Ra =
  // g*beta*(Th-Tl)*L^3/(nu*alpha),thermal expansion coefficient, thermal
  // diffusivity and kinematic viscosity T Pr;               // Prandtl number,
  // Pr = nu/alpha T Texpan_Coeff;     // 1 / K  // thermal expansion
  // coefficient T LatentHeat;       // J / g = 10^9 * mm^2 / s^2 // latent heat
  // of phase change T LatHeat_SHeatCap; // K // ratio of latent heat to
  // specific heat capacity
  // /*-----------*/
  // T Ch;           // g / mm^3 // characteristic physical high concentration
  // T Cl;           // g / mm^3 // characteristic physical low concentration
  // T CDiff;        // mm^2 / s  // concentration diffusion
  // T Cexpan_Coeff; // 1 / K // concentration expansion coefficient
  // /*----------------------------conversion
  // factors----------------------------*/
  // /*-----------*/
  // T Conv_L;       // mm
  // T Conv_Time;    // s
  // T Conv_U;       // mm / s
  // T Conv_rho;     // g / mm^3
  // T Conv_Mass;    // g
  // T Conv_VisKine; // mm^2 / s
  // T Conv_Force;   // g mm / s^2
  // T Conv_Acc;     // mm / s^2
  // T Conv_P;       // g / (mm s^2) = Pa
  // /*-----------*/
  // T Conv_dT;    // K
  // T Conv_TDiff; // mm^2 / s
  // /*-----------*/
  // T Conv_dC; // kg / m^3
  // T Conv_ConcDiff;
  // /*----------------------------lattice param----------------------------*/
  // T Lattice_charL;   // lattice domain size = Ni
  // T Lattice_charU;   // char lattice velocity
  // T Lattice_VisKine; // lattice kinematic viscosity
  // T Lattice_Re;      // grid Reynolds number
  // T Lattice_g;       // lattice gravity                 9.8 m / s^2 = 9800 mm
  // / s^2
  // /*-----------*/
  // // lattice deltaTemp = 1
  // T Lattice_gbetaT;           // lattice g * beta
  // T Lattice_betaT;            // lattice thermal expansion coefficient
  // T Lattice_TDiff;            // lattice thermal diffusion
  // T Lattice_TCond;            // lattice thermal conductivity
  // T Lattice_SHeatCap;         // lattice specific heat capacity at constant
  // pressure T Lattice_LatentHeat;       // lattice latent heat T
  // Lattice_LatHeat_SHeatCap; // lattice ratio of latent heat to specific heat
  // capacity
  // /*-----------*/
  // T Lattice_gbetaC; // lattice g * beta
  // T Lattice_betaC;  // lattice concentration expansion coefficient
  // T Lattice_CDiff;  // lattice concentration diffusion
  // /*----------------------------relaxation
  // parameters----------------------------*/ T Lattice_RT; // relaxation time
  // T OMEGA;
  // T Lattice_RTT;
  // T OMEGAT;
  // T Lattice_RTC;
  // T OMEGAC;
  // T cs2;

  virtual T GetLattice_RT() const = 0;
  virtual T GetOMEGA() const = 0;
  virtual T getLatticeRho(T rho_phys) const = 0;
  virtual T getLatRhoInit() const = 0;
  virtual T getPhysRho(T Lattice_rho) const = 0;
  virtual T GetLattice_gbeta() const { return 0; }
  virtual T getLatticeU(T U_phys) const { return 0; }
  virtual T getPhysU(T Lattice_U) const { return 0; }
};

template <typename T>
struct BaseConverter final : public AbstractConverter<T> {
  T deltaX;   // mm          // spacing between two lattice cells
  T deltaT;   // s           // time step
  T charL;    // mm          // i.e. simulation domain size
  T charU;    // mm / s      // maximal or expected velocity during simulation
  T charP;    // g / (mm s^2) = Pa
  T VisKine;  // mm^2 / s  // kinematic viscosity
  T rho;      // g / mm^3    // density
  T Re;       // Reynolds number, Re = U0*L0/nu = charU*charL/VisKine
  /*-----------*/
  T Conv_L;        // mm
  T Conv_Time;     // s
  T Conv_U;        // mm / s
  T Conv_rho;      // g / mm^3
  T Conv_Mass;     // g
  T Conv_VisKine;  // mm^2 / s
  T Conv_Force;    // g mm / s^2
  T Conv_Acc;      // mm / s^2
  T Conv_P;        // g / (mm s^2) = Pa
  /*-----------*/
  T Lattice_charL;    // lattice domain size = Ni
  T Lattice_charU;    // char lattice velocity
  T Lattice_VisKine;  // lattice kinematic viscosity
  T Lattice_Re;       // grid Reynolds number
  T Lattice_g;        // lattice gravity                 9.8 m / s^2 = 9800 mm / s^2

  T Lattice_RT;  // relaxation time
  T OMEGA;

  T cs2;

  BaseConverter(T cs2_) : cs2(cs2_) {}
  T GetLattice_RT() const override { return Lattice_RT; }
  T GetOMEGA() const override { return OMEGA; }
  T getLatticeRho(T rho_phys) const override { return rho_phys / Conv_rho; }
  T getLatRhoInit() const override { return T(1); }

  // U method
  T getLatticeU(T U_phys) const override { return U_phys / Conv_U; }
  T getPhysU(T Lattice_U) const override { return Lattice_U * Conv_U; }
  template <unsigned int D>
  Vector<T, D> getLatticeU(const Vector<T, D> &U_phys) {
    return U_phys / Conv_U;
  }
  template <unsigned int D>
  Vector<T, D> getPhysU(const Vector<T, D> &Lattice_U) {
    return Lattice_U * Conv_U;
  }

  T getPhysRho(T Lattice_rho) const override { return Lattice_rho * Conv_rho; }

  T getPhysTime(T Lattice_Time) { return Lattice_Time * deltaT; }
  /*--------------------Basic Converters--------------------*/
  void Converter(T deltaX_, T deltaT_, T rho_, T charL_, T charU_, T VisKine_, T charP_ = T(0)) {
    /*------------phys param-------------*/
    deltaX = deltaX_;
    deltaT = deltaT_;
    charL = charL_;
    charU = charU_;
    VisKine = VisKine_;
    rho = rho_;
    Re = charU_ * charL_ / VisKine_;
    /*----------conversion factors--------------*/
    Conv_L = deltaX;
    Conv_Time = deltaT;
    Conv_rho = rho;
    Conv_U = Conv_L / Conv_Time;                      // dx / dt
    Conv_Mass = Conv_rho * Conv_L * Conv_L * Conv_L;  // m = rho * dx^3
    Conv_VisKine = Conv_L * Conv_L / Conv_Time;       // nu = dx^2 / dt
    Conv_Force =
      Conv_Mass * Conv_L / Conv_Time / Conv_Time;  // F = m * dx / dt^2 = rho * dx^4 / dt^2
    Conv_Acc = Conv_L / Conv_Time / Conv_Time;     // a = dx / dt^2
    Conv_P = Conv_Force / Conv_L / Conv_L;         // P = F / dx^2 = rho * dx^2 / dt^2
    /*-----------lattice param------------*/
    Lattice_charU = charU / Conv_U;
    Lattice_charL = charL / Conv_L;
    Lattice_VisKine = VisKine / Conv_VisKine;
    Lattice_Re = Lattice_charU / Lattice_VisKine;
    Lattice_g = T(9800) / Conv_Acc;

    OMEGA = T(1) / Lattice_RT;
  }

  /*-------------------- Converters--------------------*/
  void ConvertFromRT(T deltaX_, T LatRT_, T rho_, T charL_, T charU_, T VisKine_) {
    Lattice_RT = LatRT_;
    Converter(deltaX_, (LatRT_ - T(0.5)) * cs2 * deltaX_ * deltaX_ / VisKine_, rho_, charL_, charU_,
              VisKine_);
  }
  void ConvertFromTimeStep(T deltaX_, T deltaT_, T rho_, T charL_, T charU_, T VisKine_) {
    Lattice_RT = T(0.5) + deltaT_ * VisKine_ / (cs2 * deltaX_ * deltaX_);
    Converter(deltaX_, deltaT_, rho_, charL_, charU_, VisKine_);
  }
  void SimplifiedConvertFromViscosity(int Ni_, T charU_, T VisKine_) {
    // deltaX, deltaT, rho = 1
    Lattice_RT = T(0.5) + VisKine_ / cs2;
    Converter(T(1), T(1), T(1), static_cast<T>(Ni_), charU_, VisKine_);
  }
  void SimplifiedConverterFromRT(int Ni_, T charU_, T LatRT_) {
    // deltaX, deltaT, rho = 1
    Lattice_RT = LatRT_;
    Converter(T(1), T(1), T(1), static_cast<T>(Ni_), charU_, (LatRT_ - T(0.5)) * cs2);
  }

  void check(int &check_status);
};

template <typename T>
struct TempConverter final : public AbstractConverter<T> {
  T Th;                // K      //characteristic physical high temperature
  T Tl;                // K      //characteristic physical low temperature
  T TCond;             // W / (mm K)// thermal conductivity
  T SHeatCap;          // J / (g K) = 10^9 * mm^2 / s^2 / K // specific heat capacity at
                       // constant pressure
  T TDiff;             // mm^2 / s  // thermal diffusion
  T Ra;                // Rayleigh number, Ra = g*beta*(Th-Tl)*L^3/(nu*alpha),thermal
                       // expansion coefficient, thermal diffusivity and kinematic viscosity
  T Pr;                // Prandtl number, Pr = nu/alpha
  T Texpan_Coeff;      // 1 / K  // thermal expansion coefficient
  T LatentHeat;        // J / g = 10^9 * mm^2 / s^2 // latent heat of phase change
  T LatHeat_SHeatCap;  // K // ratio of latent heat to specific heat capacity
  /*-----------*/
  T Conv_dT;     // K
  T Conv_TDiff;  // mm^2 / s
  /*-----------*/
  // lattice deltaTemp = 1
  T Lattice_gbetaT;            // lattice g * beta
  T Lattice_betaT;             // lattice thermal expansion coefficient
  T Lattice_TDiff;             // lattice thermal diffusion
  T Lattice_TCond;             // lattice thermal conductivity
  T Lattice_SHeatCap;          // lattice specific heat capacity at constant pressure
  T Lattice_LatentHeat;        // lattice latent heat
  T Lattice_LatHeat_SHeatCap;  // lattice ratio of latent heat to specific heat
                               // capacity
  T Lattice_TInit;             // lattice initial temperature
  T TInit;                     // physical initial temperature

  T Lattice_RTT;
  T OMEGAT;

  T cs2;

  BaseConverter<T> &BaseConv;
  T GetLattice_RT() const override { return Lattice_RTT; }
  T GetOMEGA() const override { return OMEGAT; }
  T GetLattice_gbeta() const override { return Lattice_gbetaT; }
  T getLatticeRho(T T_phys) const override {
    // normalized temperature
    T Lattice_T = (T_phys - Tl) / Conv_dT;
    if (Lattice_T < 0 || Lattice_T > 1) {
      std::cout << "Error: Lattice_T out of range [0,1]" << std::endl;
      exit(-1);
    }
    return Lattice_T;
  }
  T getLatRhoInit() const override { return Lattice_TInit; }
  T getPhysRho(T Lattice_T) const override { return Lattice_T * Conv_dT + Tl; }
  T getLatticeDTemp(T dT_phys) { return dT_phys / Conv_dT; }
  T getPhysDTemp(T Lattice_dT) { return Lattice_dT * Conv_dT; }
  TempConverter(T cs2_, BaseConverter<T> &baseconv, T init)
      : cs2(cs2_), BaseConv(baseconv), TInit(init) {}
  /*--------------------Converters--------------------*/
  void Converter(T Tl_, T Th_, T TDiff_) {
    Tl = Tl_;
    Th = Th_;
    TDiff = TDiff_;
    Pr = BaseConv.VisKine / TDiff_;
    /*----------conversion factors--------------*/
    Conv_dT = Th_ - Tl_;
    Conv_TDiff = BaseConv.Conv_VisKine;
    /*-----------lattice param------------*/
    Lattice_TDiff = TDiff_ / Conv_TDiff;
    Lattice_RTT = T(0.5) + BaseConv.deltaT * TDiff_ / (cs2 * BaseConv.deltaX * BaseConv.deltaX);
    OMEGAT = T(1) / Lattice_RTT;
    SetLatInit_fromPhys(TInit);
  }
  void SetLatInit_fromPhys(T physT) { Lattice_TInit = getLatticeRho(physT); }
  void SetLatInit(T latT) { Lattice_TInit = latT; }
  /*--------------------Temperature Converters--------------------*/
  bool TExpanConverted = false;
  bool LatentHeatConverted = false;
  bool SHeatCapConverted = false;
  void ConvertTempFromSHeatCap_and_TCond(T Tl_, T Th_, T TCond_, T SHeatCap_) {
    Converter(Tl_, Th_, TCond_ / SHeatCap_ / BaseConv.rho);
    TCond = TCond_;
    SHeatCap = SHeatCap_;

    SHeatCapConverted = true;
  }
  void ConvertTempFromSHeatCap_and_TCond_with_Texpan(T Tl_, T Th_, T TCond_, T SHeatCap_,
                                                     T TexpanCoeff_) {
    Converter(Tl_, Th_, TCond_ / SHeatCap_ / BaseConv.rho);
    TCond = TCond_;
    SHeatCap = SHeatCap_;
    Texpan_Coeff = TexpanCoeff_;
    Ra = T(9800) * TexpanCoeff_ * Conv_dT * pow(BaseConv.charL, 3) /
         (BaseConv.VisKine * TDiff);  // Ra = g*beta*(Th-Tl)*L^3/(nu*alpha)

    Lattice_betaT = TexpanCoeff_ * Conv_dT;
    Lattice_gbetaT = BaseConv.Lattice_g * Lattice_betaT;

    SHeatCapConverted = true;
    TExpanConverted = true;
  }
  void ConvertTempFromTDiff_with_Ra(T Tl_, T Th_, T TDiff_, T Ra_) {
    Converter(Tl_, Th_, TDiff_);
    Ra = Ra_;  // Ra = g*beta*(Th-Tl)*L^3/(nu*alpha)
    // Lattice_gbetaT = Ra * VisKine * TDiff / (charL * charL * charL * Conv_dT)
    // / Conv_Acc * Conv_dT;
    Lattice_gbetaT = Ra_ * BaseConv.VisKine * TDiff_ / pow(BaseConv.charL, 3) / BaseConv.Conv_Acc;
    Lattice_betaT = Lattice_gbetaT / BaseConv.Lattice_g;
    Texpan_Coeff = Lattice_betaT / Conv_dT;

    TExpanConverted = true;
  }
  void ConvertTempFromSHeatCap_and_TCond_with_Ra(T Tl_, T Th_, T TCond_, T SHeatCap_, T Ra_) {
    Converter(Tl_, Th_, TCond_ / SHeatCap_ / BaseConv.rho);
    TCond = TCond_;
    SHeatCap = SHeatCap_;
    Ra = Ra_;  // Ra = g*beta*(Th-Tl)*L^3/(nu*alpha)

    Lattice_gbetaT = Ra * BaseConv.VisKine * TDiff / (pow(BaseConv.charL, 3) * Conv_dT) /
                     BaseConv.Conv_Acc * Conv_dT;
    Lattice_betaT = Lattice_gbetaT / BaseConv.Lattice_g;
    Texpan_Coeff = Lattice_betaT / Conv_dT;

    SHeatCapConverted = true;
    TExpanConverted = true;
  }
  void ConvertLatentHeat(T LatentHeat_) {
    LatentHeat = LatentHeat_;  // J/g = 1e9 mm^2/s^2
    LatentHeatConverted = true;
  }
  // check
  void check(int &check_status);
};

template <typename T>
struct ConcConverter final : public AbstractConverter<T> {
  T Ch;            // g / mm^3 // characteristic physical high concentration
  T Cl;            // g / mm^3 // characteristic physical low concentration
  T CDiff;         // mm^2 / s  // concentration diffusion
  T Cexpan_Coeff;  // 1 / K // concentration expansion coefficient
  /*-----------*/
  T Conv_dC;  // kg / m^3
  T Conv_ConcDiff;
  /*-----------*/
  T Lattice_gbetaC;  // lattice g * beta
  T Lattice_betaC;   // lattice concentration expansion coefficient
  T Lattice_CDiff;   // lattice concentration diffusion
  /*-----------*/
  T Lattice_CInit;  // lattice initial concentration
  T CInit;          // physical initial concentration

  T Lattice_RTC;
  T OMEGAC;

  T cs2;
  // non-uniform time step, cause TDiff is several orders of magnitude larger
  // than CDiff ~100
  int TimeStepCoeff = 0;

  BaseConverter<T> &BaseConv;
  T GetLattice_RT() const override { return Lattice_RTC; }
  T GetOMEGA() const override { return OMEGAC; }
  T GetLattice_gbeta() const override { return Lattice_gbetaC; }
  T getLatticeDConc(T dConc_phys) { return dConc_phys / Conv_dC; }
  T getLatticeRho(T Conc_phys) const override {
    // normalized concentration
    T Lattice_Conc = (Conc_phys - Cl) / Conv_dC;
    if (Lattice_Conc < 0 || Lattice_Conc > 1) {
      std::cout << "Error: Lattice_Conc out of range [0,1]" << std::endl;
      exit(-1);
    }
    return Lattice_Conc;
  }
  T getLatRhoInit() const override { return Lattice_CInit; }
  T getPhysRho(T Lattice_Conc) const override { return Lattice_Conc * Conv_dC + Cl; }

  ConcConverter(T cs2_, BaseConverter<T> &baseconv, T init)
      : cs2(cs2_), BaseConv(baseconv), CInit(init) {}
  /*--------------------Converters--------------------*/
  void Converter(T Cl_, T Ch_, T CDiff_) {
    Cl = Cl_;
    Ch = Ch_;
    CDiff = CDiff_;
    /*----------conversion factors--------------*/
    Conv_dC = Ch_ - Cl_;
    Conv_ConcDiff = BaseConv.Conv_VisKine;
    /*-----------lattice param------------*/
    Lattice_CDiff = CDiff_ / Conv_ConcDiff;
    Lattice_RTC = T(0.5) + BaseConv.deltaT * CDiff_ / (cs2 * BaseConv.deltaX * BaseConv.deltaX);
    OMEGAC = T(1) / Lattice_RTC;
    SetLatInit_fromPhys(CInit);
  }
  void SetLatInit_fromPhys(T physC) { Lattice_CInit = getLatticeRho(physC); }
  void SetLatInit(T latC) { Lattice_CInit = latC; }

  /*--------------------Concentration Converters--------------------*/
  bool CExpanConverted = false;
  void ConvertConc(T Cl_, T Ch_, T CDiff_) { Converter(Cl_, Ch_, CDiff_); }
  void ConvertConc_withCExpan(T Cl_, T Ch_, T CDiff_, T CexpanCoeff_) {
    if (TimeStepCoeff == 0) {
      Converter(Cl_, Ch_, CDiff_);
    } else {
      Converter(Cl_, Ch_, CDiff_ * TimeStepCoeff);
    }
    Cexpan_Coeff = CexpanCoeff_;
    Lattice_betaC = Cexpan_Coeff * Conv_dC;
    Lattice_gbetaC = BaseConv.Lattice_g * Lattice_betaC;

    CExpanConverted = true;
  }
  void Enable_Non_Uniform_TimeStep(int TimeStepCoeff_ = 100) { TimeStepCoeff = TimeStepCoeff_; }
  void check(int &check_status);
};

template <typename T>
struct PhaseDiagramConverter {
  T T_Melt;      // K
  T T_Eute;      // K
  T m_Liquidus;  //
  T m_Solidus;   //
  T Part_Coef;   //

  T Lattice_T_Melt;
  T Lattice_T_Eute;
  T Lattice_m_Liq;
  T Lattice_m_Sol;

  TempConverter<T> &TempConv;
  ConcConverter<T> &ConcConv;

  PhaseDiagramConverter(TempConverter<T> &TempConv_, ConcConverter<T> &ConcConv_, T t_melt,
                        T t_eute, T m_solidus, T m_liquidus)
      : TempConv(TempConv_), ConcConv(ConcConv_), T_Melt(t_melt), T_Eute(t_eute),
        m_Solidus(m_solidus), m_Liquidus(m_liquidus) {
    Lattice_T_Melt = TempConv.getLatticeRho(T_Melt);
    Lattice_T_Eute = TempConv.getLatticeRho(T_Eute);
    Lattice_m_Sol = getLatticePhaseDiagramSlope(m_Solidus);
    Lattice_m_Liq = getLatticePhaseDiagramSlope(m_Liquidus);
    Part_Coef = Lattice_m_Liq / Lattice_m_Sol;
  }

  T getLatticePhaseDiagramSlope(T Slope_phys) {
    return Slope_phys / (TempConv.Conv_dT / ConcConv.Conv_dC);
  }
  T getTotalConcSource(T ConcLatLiq) { return ConcLatLiq * (1 - Part_Coef); }
  T getSolifiedConc(T ConcLatLiq) { return ConcLatLiq * Part_Coef; }
  // phase diagram, use normalized concentration
  T get_LatTliq(T latconc) { return Lattice_T_Melt - latconc * Lattice_m_Liq; }
  T get_LatTsol(T latconc) { return Lattice_T_Melt - latconc * Lattice_m_Sol; }
  T get_LatC_fromLiq(T latT) { return (Lattice_T_Melt - latT) / Lattice_m_Liq; }
  T get_LatC_fromSol(T latT) { return (Lattice_T_Melt - latT) / Lattice_m_Sol; }
  // phase diagram, use physical concentration
  T get_Tliq(T conc) { return T_Melt - conc * m_Liquidus; }
  T get_Tsol(T conc) { return T_Melt - conc * m_Solidus; }
  T get_C_fromLiq(T Temp) { return (T_Melt - Temp) / m_Liquidus; }
  T get_C_fromSol(T Temp) { return (T_Melt - Temp) / m_Solidus; }
  virtual void check(int &check_status) = 0;
};

template <typename T>
struct ZSConverter final : public PhaseDiagramConverter<T> {
  T GT_Coef;  // Gibbs–Thomson coefficient, mm * K

  T Lattice_GT_Coef;

  BaseConverter<T> &BaseConv;
  ZSConverter(BaseConverter<T> &BaseConv_, TempConverter<T> &TempConv_, ConcConverter<T> &ConcConv_,
              T t_melt, T t_eute, T m_solidus, T m_liquidus, T gt_coef)
      : PhaseDiagramConverter<T>(TempConv_, ConcConv_, t_melt, t_eute, m_solidus, m_liquidus),
        BaseConv(BaseConv_), GT_Coef(gt_coef) {
    Lattice_GT_Coef = GT_Coef / BaseConv.Conv_L / this->TempConv.Conv_dT;
  }
  void check(int &check_status) override;
};

template <typename T>
struct GandinConverter final : public PhaseDiagramConverter<T> {
  T DT_Mean_Bulk;   // K
  T DT_Std_Bulk;    // K
  T DT_Mean_Surf;   // K
  T DT_Std_Surf;    // K
  T Nuc_Dens_Bulk;  // mm^-3
  T Nuc_Dens_Surf;  // mm^-2
  T GrowthPara;     // mm/s/K^2

  /*-----------*/
  T Lattice_DT_Mean_Bulk;
  T Lattice_DT_Std_Bulk;
  T Lattice_DT_Mean_Surf;
  T Lattice_DT_Std_Surf;
  T Lattice_NucDens_Bulk;
  T Lattice_NucDens_Surf;
  T Lattice_GrowthPara;

  BaseConverter<T> &BaseConv;
  TempConverter<T> &TempConv;
  ConcConverter<T> &ConcConv;

  GandinConverter(BaseConverter<T> &BaseConv_, TempConverter<T> &TempConv_,
                  ConcConverter<T> &ConcConv_, T t_melt, T t_eute, T m_solidus, T m_liquidus)
      : PhaseDiagramConverter<T>(TempConv_, ConcConv_, t_melt, t_eute, m_solidus, m_liquidus),
        BaseConv(BaseConv_), TempConv(TempConv_), ConcConv(ConcConv_) {}

  void ConvertCA(T dt_mean_bulk, T dt_std_bulk, T dt_mean_surf, T dt_std_surf, T nuc_dens_bulk,
                 T nuc_dens_surf, T growthpara) {
    DT_Mean_Bulk = dt_mean_bulk;
    DT_Std_Bulk = dt_std_bulk;
    DT_Mean_Surf = dt_mean_surf;
    DT_Std_Surf = dt_std_surf;
    Nuc_Dens_Bulk = nuc_dens_bulk;
    Nuc_Dens_Surf = nuc_dens_surf;
    GrowthPara = growthpara;

    // convert to lattice UNIT
    Lattice_DT_Mean_Bulk = TempConv.getLatticeDTemp(DT_Mean_Bulk);
    Lattice_DT_Std_Bulk = TempConv.getLatticeDTemp(DT_Std_Bulk);
    Lattice_DT_Mean_Surf = TempConv.getLatticeDTemp(DT_Mean_Surf);
    Lattice_DT_Std_Surf = TempConv.getLatticeDTemp(DT_Std_Surf);
    Lattice_NucDens_Bulk = Nuc_Dens_Bulk * BaseConv.Conv_L * BaseConv.Conv_L;
    Lattice_NucDens_Surf = Nuc_Dens_Surf * BaseConv.Conv_L;
    Lattice_GrowthPara =
      GrowthPara / BaseConv.Conv_L * BaseConv.Conv_Time * TempConv.Conv_dT * TempConv.Conv_dT;
  }
  // phase diagram, use normalized concentration
  // T get_LatTliq(T conc) {
  //   return this->Lattice_T_Melt - conc * this->Lattice_m_Liq;
  // }
  // T get_LatTsol(T conc) {
  //   return this->Lattice_T_Melt - conc * this->Lattice_m_Sol;
  // }
  T EuteCorrection(T temp) {
    if (temp < this->Lattice_T_Eute)
      return this->Lattice_T_Eute;
    else
      return temp;
  }
  void check(int &check_status) override;
};

template <typename T>
class UnitConvManager {
  /*Unit converter for LB*/
  // Unit phys = Unit LB * Conversionfactor
  // pressure and temperature: shift the physical values by a characteristic
  // value -> lattice p and T ranging from 0 to 1 e.g. physT - charPhysT =
  // Conversionfactor * LBT i.e. T - T0 = (Th - T0) * Theta
 public:
  BaseConverter<T> *BaseConv;
  TempConverter<T> *TempConv;
  ConcConverter<T> *ConcConv;
  PhaseDiagramConverter<T> *CAConv;
  std::vector<AbstractConverter<T> *> ConvList;
  UnitConvManager(AbstractConverter<T> *convlist) : ConvList(convlist) {}
  UnitConvManager(BaseConverter<T> *BaseConv_, TempConverter<T> *TempConv_ = nullptr,
                  ConcConverter<T> *ConcConv_ = nullptr,
                  PhaseDiagramConverter<T> *CAConv_ = nullptr)
      : BaseConv(BaseConv_), TempConv(TempConv_), ConcConv(ConcConv_), CAConv(CAConv_) {
    ConvList.push_back(BaseConv);
    if (TempConv != nullptr) {
      ConvList.push_back(TempConv);
    }
    if (ConcConv != nullptr) {
      ConvList.push_back(ConcConv);
    }
  }
  void Check_and_Print();
};

// refine converter
template <typename T>
struct RefineConverter {
  // physical relaxation time(RT) reamains unchanged between different refinement levels
  // (Lat_RTC - 0.5) * deltaTC = (LatRTF - 0.5) * deltaTF
  // Lat_RTC = 0.5*Lat_RTF + 0.25
  // Lat_RTF = 2 * Lat_RTC - 0.5
  // OmegaC = 1 / (0.5 / OmegaF + 0.25)
  // OmegaF = 1 / (2 / OmegaC - 0.5)
  static inline T getOmegaF(T omegaC) { return T(1) / (T(2) / omegaC - T(0.5)); }
  static inline T getOmegaC(T omegaF) { return T(1) / (T(0.5) / omegaF + T(0.25)); }
  // recursively get OmegaF
  static T getOmegaF(T omegaC, std::uint8_t level) {
    if (level == std::uint8_t(0)) return omegaC;
    if (level == std::uint8_t(1)) return getOmegaF(omegaC);
    return getOmegaF(getOmegaF(omegaC, level - 1));
  }
  // recursively get OmegaC
  static T getOmegaC(T omegaF, std::uint8_t level) {
    if (level == std::uint8_t(0)) return omegaF;
    if (level == std::uint8_t(1)) return getOmegaC(omegaF);
    return getOmegaC(getOmegaC(omegaF, level - 1));
  }

  // continuous distribution function remains unchanged between different refinement levels
  // gC + (0.5*1/Lat_RTC)(geq - gC) = gF + (0.5*1/Lat_RTF)(geq - gF)
  // gC = gF + (0.5*1/Lat_RTF)(geq - gF) = gF + (0.5*OmegaF)(gF - geq)
  // gF = gC - (0.25*1/Lat_RTC)(geq - gC) = gC - (0.25*OmegaC)(gC - geq)
  static inline T getPopF(T popC, T popeq, T omegaC) {
    return popC - T(0.25) * omegaC * (popC - popeq);
  }
  static inline T getPopC(T popF, T popeq, T omegaF) {
    return popF + T(0.5) * omegaF * (popF - popeq);
  }
};
