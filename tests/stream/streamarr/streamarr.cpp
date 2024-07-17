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
int shift;
int rotatecount;

void readParam() {
  iniReader param_reader("streamarr.ini");
  N = param_reader.getValue<int>("Mesh", "N");
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

void set(StreamArray<T> &arr) {
  for(int i = 0; i < arr.size(); i++) {
    arr[i] = i+1;
  }
}

void set(CyclicArray<T> &arr) {
  for(int i = 0; i < arr.size(); i++) {
    arr[i] = i+1;
  }
}

int main() {

  
  readParam();

  StreamArray<T> sarr(N, T{});
  CyclicArray<T> carr(N, T{});

  int offset = int(std::sqrt(N));

  sarr.copyToBack();
  sarr.setOffset(offset);

  // std::cout << "----------------" << std::endl;

  Timer MainLoopTimer;
  MainLoopTimer.START_TIMER();
  for(int i = 0; i < rotatecount; i++) {
    set(sarr);
    sarr.rotate();
  }
  MainLoopTimer.END_TIMER();
  std::cout << "StreamArray efficiency: " << MainLoopTimer.GetDurationCount_Only() <<  std::endl;

  MainLoopTimer.START_TIMER();
  for(int i = 0; i < rotatecount; i++) {
    set(carr);
    carr.rotate(offset);
  }
  MainLoopTimer.END_TIMER();
  std::cout << "CyclicArray efficiency: " << MainLoopTimer.GetDurationCount_Only() <<  std::endl;

  

  return 0;
}