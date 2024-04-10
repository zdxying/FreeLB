// a test file to find the most efficient way to get Remain sites
//

#include "backerase.h"
#include "backerase.hh"
#include "freelb.h"
#include "freelb.hh"
#include "geometry_.h"
#include "geometry_.hh"

using T = FLOAT;

int Ni;
int Nj;
int MaxStep;
int OutputStep;

void readParam() {
  /*reader*/
  iniReader param_reader("test.ini");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  MaxStep = param_reader.getValue<int>("Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Settings", "OutputStep");
  // std::cout << "[Settings]:"
  //           << "total N:         " << Ni * Nj << "\n"
  //           << "TotalStep:         " << MaxStep << "\n"
  //           << "OutputStep:        " << OutputStep << "\n"
  //           << "----------------------------------------------\n"
  //           << std::endl;
}

int main() {
  readParam();
  Geometry2D_legacy<T> Geo(Ni, Nj);
  Geo.circle(Ni / 10, Ni / 2, Nj / 2, 1);

  Test_Thread Test(Ni, Nj, Geo.Flag);
  Timer timer1;
  Test.Get_Remain_pushback();
  std::cout << "Get_Remain_pushback: " << timer1.GetTimeElapsed() << std::endl;
  std::cout << "RemainSize: " << Test.Remain.size() << std::endl;
  timer1.START_TIMER();
  Test.expand_Get_Remain_pushback();
  std::cout << "expand_Get_Remain_pushback: " << timer1.GetTimeElapsed()
            << std::endl;
  std::cout << "RemainSize: " << Test.Remain.size() << std::endl;
  
  Test.getremain<&Test_Thread::get>();

  // Test.Get_Remain_erase();
  // std::cout << "Get_Remain_erase: " << timer1.GetTimeElapsed() << std::endl;
  // std::cout << "RemainSize: " << Test.Remain.size() << std::endl;

  std::cout << "[details]: " << std::endl;
  std::cout << "pushback_time: " << Test.pushback_time << std::endl;
  std::cout << "pushback_merge_time: " << Test.pushback_merge_time << std::endl;
  std::cout << "expand_pushback_time: " << Test.expand_pushback_time
            << std::endl;
  std::cout << "expand_pushback_merge_time: " << Test.expand_pushback_merge_time
            << std::endl;
  // std::cout << "erase_time: " << Test.erase_time << std::endl;
  // std::cout << "erase_merge_time: " << Test.erase_merge_time << std::endl;
  // std::cout << "erase_divide_time: " << Test.erase_divide_time << std::endl;

  // CHECK
  // std::cout << "Remain_Count: " << Test.Remain_Count << std::endl;
  // std::cout << "Erase_Count:  " << Test.Erase_Count << std::endl;
  // test vectors
  // remain
}