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

#include "freelb.h"
#include "freelb.hh"

#include "template.h"

#define TEMPGEN_LATSETS tempgen::D2Q5, tempgen::D2Q9, tempgen::D3Q7, tempgen::D3Q15, tempgen::D3Q19, tempgen::D3Q27

int main() {
  std::string outdir = "./output/";
  std::string momentfname = outdir + "cse.h";
  std::string feqfname = outdir + "feq.h";
  std::string forcefname = outdir + "force.h";

  DirCreator::Create_Dir(outdir);

  tempgen::rhoImplgen<TEMPGEN_LATSETS> rho(momentfname);
  tempgen::sourcerhoImplgen<TEMPGEN_LATSETS> sourcerho(momentfname);
  tempgen::UImplgen<TEMPGEN_LATSETS> U(momentfname);
  tempgen::forceUImplgen<TEMPGEN_LATSETS> forceU(momentfname);
  tempgen::rhoUImplgen<TEMPGEN_LATSETS> rhoU(momentfname);
  tempgen::forcerhoUImplgen<TEMPGEN_LATSETS> forcerhoU(momentfname);
  tempgen::Pi_ab_neqgen<TEMPGEN_LATSETS> Pi_ab_neq(momentfname);
  tempgen::forcePi_ab_neqgen<TEMPGEN_LATSETS> forcePi_ab_neq(momentfname);
  tempgen::stressgen<TEMPGEN_LATSETS> stress(momentfname);
  tempgen::strainRategen<TEMPGEN_LATSETS> strainRate(momentfname);
  tempgen::shearRateMagImplgen<TEMPGEN_LATSETS> shearRateMag(momentfname);

  tempgen::SecondOrderImplgen<TEMPGEN_LATSETS> SecondOrder(feqfname);

  tempgen::ForcePopImplgen<TEMPGEN_LATSETS> ForcePop(forcefname);
  tempgen::ScalarForcePopImplgen<TEMPGEN_LATSETS> ScalarForcePop(forcefname);

  rho.generateAll();
  sourcerho.generateAll();
  U.generateAll();
  forceU.generateAll();
  rhoU.generateAll();
  forcerhoU.generateAll();
  Pi_ab_neq.generateAll();
  forcePi_ab_neq.generateAll();
  stress.generateAll();
  strainRate.generateAll();
  shearRateMag.generateAll();

  SecondOrder.generateAll();

  ForcePop.generateAll();
  ScalarForcePop.generateAll();


}
