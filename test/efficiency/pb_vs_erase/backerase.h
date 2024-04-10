// header file for test.cpp
#pragma once

#include <omp.h>

#include <iostream>
#include <vector>

#include "utils/timer.h"
#include "utils/util.h"

class Test_Thread {
 public:
  int Ni;
  int Nj;
  int N;
  int *Flag;

  std::vector<int> Remain;
#ifdef _REMAIN_PUSHBACK_P
  std::vector<std::vector<int>> Remain_Thread;
#endif
  int Erase_Count = 0;

  // constructor
  Test_Thread(int ni, int nj, int *flag);

  void Get_Remain_pushback();
  double pushback_time = 0;
  double pushback_merge_time = 0;
  void expand_Get_Remain_pushback();
  double expand_pushback_time = 0;
  double expand_pushback_merge_time = 0;
  void Get_Remain_erase();
  double erase_time = 0;
  double erase_merge_time = 0;
  double erase_divide_time = 0;

  template <void (Test_Thread::*get)(int)>
  void getremain() {
    Index2D::Traverse_Bulk1<get>(Ni, Nj, 0);
  }

//
#ifdef _REMAIN_PUSHBACK_P
  void get_p(int thread_num, int id) {
    if (Flag[id] == 0) {
      Remain_Thread[thread_num].emplace_back(id);
    }
  }
#else
  void get(int id) {
    if (Flag[id] == 0) {
      Remain.emplace_back(id);
    }
  }
#endif

  // function object
  //   struct Get_Remain {
  //     #ifdef _REMAIN_PUSHBACK_P
  //     void operator()(int thread_num, int id) const {
  //     if (Flag[id] == 0) {
  //       Remain_Thread[thread_num].emplace_back(id);
  //     }
  //   }
  //   #else
  //   void operator()(int id) const {
  //     if (Flag[id] == 0) {
  //       Remain.emplace_back(id);
  //     }
  //   }

  //   #endif
  // };
};