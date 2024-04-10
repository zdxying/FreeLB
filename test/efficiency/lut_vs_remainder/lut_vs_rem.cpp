// test files for the LUT vs remainder
// lookuptable vs remainder

#include <vector>

#include "data_struct/Vector.h"
#include "utils/timer.h"

template <unsigned int D>
class Geometry {
 private:
  int _Nx;
  int _Ny;
  int _Nz;
  int N;
  int thread_num = 4;
  int NxNy;
  // look up table
  std::vector<Vector<int, D>> _LUT;
  std::vector<int> Field;

 public:
  Geometry(int Nx, int Ny, int Nz = 1) : _Nx(Nx), _Ny(Ny), _Nz(Nz) {
    N = _Nx * _Ny * _Nz;
    NxNy = _Nx * _Ny;
    _LUT.reserve(N);
    Field.resize(N, 0);
    InitLut();
  }

  void InitLut() {
    if constexpr (D == 3) {
      for (int k = 0; k < _Nz; ++k) {
        for (int j = 0; j < _Ny; ++j) {
          for (int i = 0; i < _Nx; ++i) {
            _LUT.emplace_back(i, j, k);
          }
        }
      }
    } else if constexpr (D == 2) {
      for (int j = 0; j < _Ny; ++j) {
        for (int i = 0; i < _Nx; ++i) {
          _LUT.emplace_back(i, j);
        }
      }
    }
  }

  Vector<int, D> RetLoc(int id) const {
    Vector<int, D> loc;
    if constexpr (D == 2) {
      loc[0] = id % _Nx;
      loc[1] = id / _Nx;
    } else if constexpr (D == 3) {
      int temp = id;
      loc[0] = temp % _Nx;
      temp /= _Nx;
      loc[1] = temp % _Ny;
      loc[2] = temp / _Ny;
    }
    return loc;
  }

  void GetLoc(int id, Vector<int, D>& loc) const {
    if constexpr (D == 2) {
      loc[0] = id % _Nx;
      loc[1] = id / _Nx;
    } else if constexpr (D == 3) {
      int temp = id;
      loc[0] = temp % _Nx;
      temp /= _Nx;
      loc[1] = temp % _Ny;
      loc[2] = temp / _Ny;
    }
  }

  Vector<int, D> RetLoc_(int id) const {
    Vector<int, D> loc;
    GetLoc_(id, loc);
    return loc;
  }
  
  void GetLoc_(int id, Vector<int, D>& loc) const {
    if constexpr (D == 2) {
      loc[1] = id / _Nx;
      loc[0] = id - loc[1] * _Nx;
    } else if constexpr (D == 3) {
      loc[2] = id / NxNy;
      int temp = id - loc[2] * NxNy;
      loc[1] = temp / _Nx;
      loc[0] = temp - loc[1] * _Nx;
    }
  }

  void Lut(int step) {
    int i = 0;
    while (i < step) {
      #pragma omp parallel for num_threads(thread_num) schedule(static)
      for (int id = 0; id < N; ++id) {
        Vector<int, D> loc = _LUT[id];
        int id_ = loc[0] + loc[1] * _Nx;
        Field[id_] = id;
      }
      ++i;
    }
  }

  void Rem(int step) {
    int i = 0;
    while (i < step) {
      #pragma omp parallel for num_threads(thread_num) schedule(static)
      for (int id = 0; id < N; ++id) {
        Vector<int, D> loc;
        GetLoc(id, loc);
        int id_ = loc[0] + loc[1] * _Nx;
        Field[id_] = id;
      }
      ++i;
    }
  }

  void Rem_(int step) {
    int i = 0;
    while (i < step) {
      #pragma omp parallel for num_threads(thread_num) schedule(static)
      for (int id = 0; id < N; ++id) {
        // Vector<int, D> loc;
        // GetLoc_(id, loc);
        const Vector<int, D>& loc = RetLoc_(id);
        int id_ = loc[0] + loc[1] * _Nx;
        Field[id_] = id;
      }
      ++i;
    }
  }
};


int main() {
  const int Nx = 100;
  const int Ny = 100;
  const int Nz = 100;

  int step = 100;

  Geometry<2> geo2d(Nx, Ny);
  // std::vector<Vector<int, 2>> testvec2d;
  // std::vector<Vector<int, 2>> testvec2d_;
  Geometry<3> geo3d(Nx, Ny, Nz);
  // std::vector<Vector<int, 3>> testvec3d;
  // std::vector<Vector<int, 3>> testvec3d_;

  // testvec2d.reserve(Nx * Ny);
  // testvec3d.reserve(Nx * Ny * Nz);
  // testvec2d_.reserve(Nx * Ny);
  // testvec3d_.reserve(Nx * Ny * Nz);

  Timer timer;

  geo2d.Lut(step * Nz);
  double lut2d = timer.GetTimeElapsed();

  timer.START_TIMER();
  geo2d.Rem(step * Nz);
  double rem2d = timer.GetTimeElapsed();

  timer.START_TIMER();
  geo2d.Rem_(step * Nz);
  double rem2d_ = timer.GetTimeElapsed();

  timer.START_TIMER();
  geo3d.Lut(step);
  double lut3d = timer.GetTimeElapsed();

  timer.START_TIMER();
  geo3d.Rem(step);
  double rem3d = timer.GetTimeElapsed();

  timer.START_TIMER();
  geo3d.Rem_(step);
  double rem3d_ = timer.GetTimeElapsed();

  std::cout << "2D LUT: " << lut2d << " s"
            << "\n"
            << "2D REM: " << rem2d << " s"
            << "\n"
            << "2D REM_ " << rem2d_ << " s"
            << "\n"
            << "3D LUT: " << lut3d << " s"
            << "\n"
            << "3D REM: " << rem3d << " s"
            << "\n"
            << "3D REM_ " << rem3d_ << " s" << std::endl;


  return 0;
}