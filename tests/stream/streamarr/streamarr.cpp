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

using T = FLOAT;


/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
int N;
int Ni;
int Nj;
int Nk;
int shift;
int rotatecount;

void readParam() {
  iniReader param_reader("streamarr.ini");
  N = param_reader.getValue<int>("Mesh", "N");
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Nk = param_reader.getValue<int>("Mesh", "Nk");
  shift = param_reader.getValue<int>("Mesh", "shift");
  rotatecount = param_reader.getValue<int>("Mesh", "rotatecount");
}

// void print(StreamArray<T> &arr) {
//   for(int i = 0; i < arr.size(); i++) {
//     std::cout << arr[i] << " ";
//   }
//   std::cout << std::endl;
// }

// void printprev(StreamArray<T> &arr) {
//   for(int i = 0; i < arr.size(); i++) {
//     std::cout << arr.getPrevious(i) << " ";
//   }
//   std::cout << std::endl;
// }

template <typename T, unsigned int q>
using PopStreamField = GenericArrayField<StreamArray<T>, q>;
template <typename T, unsigned int q>
using PopCyclicField = GenericArrayField<CyclicArray<T>, q>;

void set(StreamArray<T>& arr) {
  for (std::size_t i = 0; i < arr.size(); i++) {
    arr[i] = i + 1;
  }
}
void set(CyclicArray<T>& arr) {
  for (std::size_t i = 0; i < arr.size(); i++) {
    arr[i] = i + 1;
  }
}

template <typename T, unsigned int q>
void Stream(PopStreamField<T, q>& Pop, const std::array<int, q>& Delta_Index) {
  for (unsigned int k = 0; k < q; ++k) {
    Pop.getField(k).rotate(Delta_Index[k]);
  }
}
template <typename T, unsigned int q>
void Stream(PopCyclicField<T, q>& Pop, const std::array<int, q>& Delta_Index) {
  for (unsigned int k = 0; k < q; ++k) {
    Pop.getField(k).rotate(Delta_Index[k]);
  }
}

template <typename T, unsigned int q>
void set(PopStreamField<T, q>& Pop) {
  for (unsigned int k = 0; k < q; ++k) {
    set(Pop.getField(k));
  }
}
template <typename T, unsigned int q>
void set(PopCyclicField<T, q>& Pop) {
  for (unsigned int k = 0; k < q; ++k) {
    set(Pop.getField(k));
  }
}


using LatSet0 = D2Q9<T>;
using LatSet1 = D3Q19<T>;

int main() {
  readParam();

  std::size_t N0 = 100 * Ni * Nj;
  std::size_t N1 = Ni * Nj * Nk;

  Vector<int, LatSet0::d> Proj0{1, 10 * Ni};
  Vector<int, LatSet1::d> Proj1{1, Ni, Ni * Nj};

  std::array<int, LatSet0::q> Delta_Index0 =
    make_Array<int, LatSet0::q>([&](int i) { return LatSet0::c[i] * Proj0; });
  std::array<int, LatSet1::q> Delta_Index1 =
    make_Array<int, LatSet1::q>([&](int i) { return LatSet1::c[i] * Proj1; });


  PopStreamField<T, LatSet0::q> PopS0(N0, T{});
  PopCyclicField<T, LatSet0::q> PopC0(N0, T{});
  PopStreamField<T, LatSet1::q> PopS1(N1, T{});
  PopCyclicField<T, LatSet1::q> PopC1(N1, T{});

  // std::cout << "----------------" << std::endl;

  Timer MainLoopTimer;
  MainLoopTimer.START_TIMER();
  for (int i = 0; i < rotatecount; i++) {
    set(PopS0);
    Stream<T, LatSet0::q>(PopS0, Delta_Index0);
  }
  MainLoopTimer.END_TIMER();
  std::cout << "StreamField2D efficiency: " << MainLoopTimer.GetDurationCount_Only()
            << std::endl;

  MainLoopTimer.START_TIMER();
  for (int i = 0; i < rotatecount; i++) {
    set(PopC0);
    Stream<T, LatSet0::q>(PopC0, Delta_Index0);
  }
  MainLoopTimer.END_TIMER();
  std::cout << "CyclicField2D efficiency: " << MainLoopTimer.GetDurationCount_Only()
            << std::endl;

  MainLoopTimer.START_TIMER();
  for (int i = 0; i < rotatecount; i++) {
    set(PopS1);
    Stream<T, LatSet1::q>(PopS1, Delta_Index1);
  }
  MainLoopTimer.END_TIMER();
  std::cout << "StreamField3D efficiency: " << MainLoopTimer.GetDurationCount_Only()
            << std::endl;

  MainLoopTimer.START_TIMER();
  for (int i = 0; i < rotatecount; i++) {
    set(PopC1);
    Stream<T, LatSet1::q>(PopC1, Delta_Index1);
  }
  MainLoopTimer.END_TIMER();
  std::cout << "CyclicField3D efficiency: " << MainLoopTimer.GetDurationCount_Only()
            << std::endl;

  return 0;
}