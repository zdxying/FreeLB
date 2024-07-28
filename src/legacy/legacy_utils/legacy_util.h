// legacy_util.h

#pragma once

#include <vector>
#include <omp.h>
#include <math.h>

template <typename T>
struct Vect2D {
  static T dot(const T *a, const int *b) { return a[0] * b[0] + a[1] * b[1]; }
  static T dot(const T *a, const T *b) { return a[0] * b[0] + a[1] * b[1]; }
  // to be compatible with old code
  static T dot(const Vector<int, 2> &a, const T *b) {
    return a[0] * b[0] + a[1] * b[1];
  }
  static T dot(const T *a, const Vector<int, 2> &b) {
    return a[0] * b[0] + a[1] * b[1];
  }
  //
  static T cross(const T *a, const int *b) { return a[0] * b[1] - b[0] * a[1]; }
  static T cross(const T *a, const T *b) { return a[0] * b[1] - b[0] * a[1]; }
  static T sqr(const T *a) { return a[0] * a[0] + a[1] * a[1]; }
  static T norm(const T *a) { return sqrt(sqr(a)); }
  static int quad(const T *a) {
    // 2 | 1
    // -----
    // 3 | 4
    if (a[0] >= 0 && a[1] >= 0)
      return 0;
    else if (a[0] < 0 && a[1] >= 0)
      return 1;
    else if (a[0] < 0 && a[1] < 0)
      return 2;
    else
      return 3;
  }
  // get relative location of point (x, y) to point (x0, y0) in a rotated
  // coordinate system，theta is the angle of rotation counterclockwise
  // call: R_Loc(x, y(cell centre), x0, y0(growth centre), theta(orinetation),
  // loc)
  static void R_Loc(T x, T y, T x0, T y0, T theta, T *loc) {
    T dx = x - x0;
    T dy = y - y0;
    loc[0] = dx * cos(theta) + dy * sin(theta);
    loc[1] = -dx * sin(theta) + dy * cos(theta);
  }
};

// index
struct Index2D {
  // return i + j * Ni
  static inline int GetId(int i, int j, int Ni) { return j * Ni + i; }

  template <typename T, typename Func>
  static inline void Traverse_Vector(const std::vector<T> &vec, const Func& func) {
    for (int i = 0; i < vec.size(); i++) {
      func(vec[i]);
    }
  }
  // traverse bulk cells with an offset, takes a lambda funcrion, serilized
  // how to call: Traverse_Bulk(Ni, Nj, offset, [](int id){func(id);});
  // note that explicitly specify the template parameter Func is not needed
  // Attention:
  // if func is very small, and called frequently, it's recommended to inline
  // the code of func to the loop body, not to use this function
  template <typename Func>
  static void Traverse_Bulk(int Ni, int Nj, int offset, const Func& func) {
    int id;
    if (offset == 0) {
      int N = Ni * Nj;
      for (int id = 0; id < N; id++) {
        func(id);
      }
    } else {
      for (int j = offset; j < Nj - offset; j++) {
        for (int i = offset; i < Ni - offset; i++) {
          id = GetId(i, j, Ni);
          func(id);
        }
      }
    }
  }

  // traverse bulk cells with an offset, takes a lambda function, omp enabled
  // how to call: Traverse_Bulk_OMP(Ni, Nj, offset, [](int id){func(id);});
  // note that explicitly specify the template parameter Func is not needed
  // Attention:
  // 1. if omp_get_thread_num() is used in func, it's recommended to use
  // Traverse_Chunk_OMP instead
  // 2. if func is small, it's recommended to use Traverse_Bulk instead
  template <typename Func>
  static void Traverse_Bulk_OMP(int Ni, int Nj, int offset, const Func& func) {
    int id;
#pragma omp parallel for private(id) num_threads(Thread_Num)
    for (int j = offset; j < Nj - offset; j++) {
      for (int i = offset; i < Ni - offset; i++) {
        id = GetId(i, j, Ni);
        func(id);
      }
    }
  }

  // traverse bulk cells with an offset, takes a lambda function, omp enabled
  // partition the bulk cells into Thread_Num chunks
  // Traverse_Chunk_OMP(Ni, Nj, offset, [](int thn, int id){func(thn,id);});
  // note that explicitly specify the template parameter Func is not needed
  //   template <typename Func>
  //   static void Traverse_Chunk_OMP(int Ni, int Nj, int offset, const Func& func) {
  //     int id, start, end, th_num, i, j;
  //     int chunk = (Nj - 2 * offset) / Thread_Num;
  // #pragma omp parallel private(id, start, end, th_num, i, j) \
//     num_threads(Thread_Num)
  //     {
  //       th_num = omp_get_thread_num();
  //       start = offset + th_num * chunk;
  //       end = (th_num == Thread_Num - 1) ? Nj - offset : start + chunk;
  //       for (j = start; j < end; j++) {
  //         for (i = offset; i < Ni - offset; i++) {
  //           id = GetId(i, j, Ni);
  //           func(th_num, id);
  //         }
  //       }
  //     }
  //   }

  // traverse peripheral/boundary cells, takes a lambda funcrion, serilized
  // how to call: Traverse_Peripheral(Ni, Nj, [](int id){func(id);});
  // offset = 0: peripheral cells; offset = 1: boundary cells
  // note that explicitly specify the template parameter Func is not needed
  // Attention:
  // if func is very small, it's recommended to inline the code of func to the
  // loop body, not to use this function
  template <typename Func>
  static void Traverse_Peripheral(int Ni, int Nj, int offset, const Func& func) {
    int id;
    // top and bottom boundary j = 0 + offset, Nj - 1 - offset
    for (int i = offset; i < Ni - offset; i++) {
      id = GetId(i, offset, Ni);
      func(id);
    }
    for (int i = offset; i < Ni - offset; i++) {
      id = GetId(i, Nj - 1 - offset, Ni);
      func(id);
    }
    // left and right boundary i = 0 + offset, Ni - 1 - offset
    for (int j = 1 + offset; j < Nj - 1 - offset; j++) {
      id = GetId(offset, j, Ni);
      func(id);
    }
    for (int j = 1 + offset; j < Nj - 1 - offset; j++) {
      id = GetId(Ni - 1 - offset, j, Ni);
      func(id);
    }
  }
  // traverse peripheral/boundary cells, takes a lambda funcrion, omp enabled
  // how to call: Traverse_Peripheral_OMP(Ni, Nj, [](int id){func(id);});
  // offset = 0: peripheral cells; offset = 1: boundary cells
  // note that explicitly specify the template parameter Func is not needed
  // Attention: omp enabled version may be less efficient than serilized
  // version
  template <typename Func>
  static void Traverse_Peripheral_OMP(int Ni, int Nj, int offset, const Func& func) {
    int id;
    // top and bottom boundary j = 0 + offset, Nj - 1 - offset
#pragma omp parallel for private(id) num_threads(Thread_Num)
    for (int i = offset; i < Ni - offset; i++) {
      id = GetId(i, offset, Ni);
      func(id);
    }
#pragma omp parallel for private(id) num_threads(Thread_Num)
    for (int i = offset; i < Ni - offset; i++) {
      id = GetId(i, Nj - 1 - offset, Ni);
      func(id);
    }
    // left and right boundary i = 0 + offset, Ni - 1 - offset
#pragma omp parallel for private(id) num_threads(Thread_Num)
    for (int j = 1 + offset; j < Nj - 1 - offset; j++) {
      id = GetId(offset, j, Ni);
      func(id);
    }
#pragma omp parallel for private(id) num_threads(Thread_Num)
    for (int j = 1 + offset; j < Nj - 1 - offset; j++) {
      id = GetId(Ni - 1 - offset, j, Ni);
      func(id);
    }
  }
};
