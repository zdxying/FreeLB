// backerase.hh

#include "backerase.h"

Test_Thread::Test_Thread(int ni, int nj, int *flag)
    : Ni(ni), Nj(nj), N(ni * nj), Flag(flag) {
  Remain.reserve(N);
#ifdef _REMAIN_PUSHBACK_P
  OMP_helper::Setup_Thread_Vec<int>(Remain_Thread, N);
  std::cout << "size of Remain_Thread: " << Remain_Thread.size() << std::endl;
#endif
}

void Test_Thread::Get_Remain_pushback() {
  Remain.clear();
  Timer timer1;
#ifdef _REMAIN_PUSHBACK_P
  // Index2D::Traverse_Chunk_OMP(Ni, Nj, 0, [this](int thread_num, int id) {
  //   if (Flag[id] == 0) {
  //     Remain_Thread[thread_num].emplace_back(id);
  //     //   remain_count[omp_get_thread_num()]++;
  //   }
  // });

#else
  Index2D::Traverse_Bulk(Ni, Nj, 0, [this](int id) {
    if (Flag[id] == 0) {
      Remain.emplace_back(id);
    }
  });
  // Index2D::Traverse_Bulk(Ni, Nj, 0, get);

#endif

  pushback_time = timer1.GetTimeElapsed();

#ifdef _REMAIN_PUSHBACK_P
  std::cout << "using parallel pushback" << std::endl;
  timer1.START_TIMER();
  OMP_helper::Merge_and_Clear<int>(Remain, Remain_Thread);
  pushback_merge_time = timer1.GetTimeElapsed();
//   for (int i = 0; i < Thread_Num; i++) {
//     Remain_Count += remain_count[i];
//   }
#endif
}

inline void Test_Thread::expand_Get_Remain_pushback() {
  Remain.clear();
  int offset = 0;
  Timer timer1;
#ifdef _REMAIN_PUSHBACK_P
  int id, start, end, th_num, i, j;
  int chunk = (Nj - 2 * offset) / Thread_Num;
#pragma omp parallel private(id, start, end, th_num, i, j) \
    num_threads(Thread_Num)
  {
    th_num = omp_get_thread_num();
    start = offset + th_num * chunk;
    end = (th_num == Thread_Num - 1) ? Nj - offset : start + chunk;
    for (j = start; j < end; j++) {
      for (i = offset; i < Ni - offset; i++) {
        id = Index2D::GetId(i, j, Ni);
        if (Flag[id] == 0) {
          Remain_Thread[th_num].emplace_back(id);
          //   remain_count[omp_get_thread_num()]++;
        }
      }
    }
  }
#else
  int id;
  for (int j = offset; j < Nj - offset; j++) {
    for (int i = offset; i < Ni - offset; i++) {
      id = Index2D::GetId(i, j, Ni);
      if (Flag[id] == 0) {
        Remain.emplace_back(id);
        // Remain_Count++;
      }
    }
  }
#endif

  expand_pushback_time = timer1.GetTimeElapsed();
#ifdef _REMAIN_PUSHBACK_P
  std::cout << "using parallel pushback" << std::endl;
  timer1.START_TIMER();
  OMP_helper::Merge_and_Clear<int>(Remain, Remain_Thread);
  expand_pushback_merge_time = timer1.GetTimeElapsed();
#endif
}

void Test_Thread::Get_Remain_erase() {
  int Id = 0;
#ifdef _REMAIN_ERASE_P
  //   int erase_count_[Thread_Num] = {0};
  std::cout << "using parallel erase" << std::endl;
  Timer timer1;
  Timer timer2;
  std::vector<std::vector<int>> Remain_Thread_;
  OMP_helper::Divide<int>(Remain, Remain_Thread_);
  erase_divide_time = timer2.GetTimeElapsed();
  std::vector<int>::iterator iter = Remain_Thread_[0].begin();
#pragma omp parallel for private(Id, iter) num_threads(Thread_Num)
  for (int i = 0; i < Remain_Thread_.size(); i++) {
    iter = Remain_Thread_[i].begin();
    while (iter != Remain_Thread_[i].end()) {
      Id = *iter;
      if (Flag[Id] == 0) {
        iter = Remain_Thread_[i].erase(iter);
        // erase_count_[omp_get_thread_num()]++;
      } else {
        iter++;
      }
    }
  }
  erase_time = timer1.GetTimeElapsed();
  timer1.START_TIMER();
  Remain.clear();
  OMP_helper::Merge_and_Clear<int>(Remain, Remain_Thread_);
  erase_merge_time = timer1.GetTimeElapsed();
//   for (int i = 0; i < Thread_Num; i++) {
//     Erase_Count += erase_count_[i];
//   }
#else
  Timer timer1;
  std::vector<int>::iterator iter = Remain.begin();
  while (iter != Remain.end()) {
    Id = *iter;
    if (Flag[Id] == 0) {
      iter = Remain.erase(iter);
      //   Erase_Count++;
    } else {
      iter++;
    }
  }
  erase_time = timer1.GetTimeElapsed();
#endif
}
