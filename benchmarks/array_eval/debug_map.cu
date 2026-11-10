// debug: dump the StreamMapArray device buffer layout to find where the
// wrapped view reads zeros from
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

using T = float;

#define CUDA_CHECK(call)                                                        \
  do {                                                                          \
    const cudaError_t err = (call);                                             \
    if (err != cudaSuccess) {                                                   \
      std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(err), __FILE__, \
                  __LINE__);                                                    \
      std::exit(1);                                                             \
    }                                                                           \
  } while (0)

__global__ void rotateLoopKernel(cudev::StreamMapArray<T>* a, int steps) {
  if (blockIdx.x != 0 || threadIdx.x != 0) return;
  for (int s = 0; s < steps; ++s) a->rotate();
}

__global__ void dumpViewKernel(const cudev::StreamMapArray<T>* a, std::size_t n,
                               T* out) {
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) out[i] = (*a)[i];
}

// host-side POD with the same layout as cudev::StreamMapArray<T>
// (count*, _data, _shift*, _start, Offset*) -- memcpy only, no ctor calls
struct MirrorPOD {
  std::size_t* count;
  T* _data;
  std::ptrdiff_t* _shift;
  T* _start;
  std::ptrdiff_t* Offset;
};

int main() {
  const std::size_t n = 1000;
  StreamMapArray<T> a(n, T{});
  for (std::size_t i = 0; i < n; ++i) a[i] = T(i) + T(1);  // content = i+1
  a.setOffset(3);
  a.copyToDevice();

  MirrorPOD m;
  static_assert(sizeof(MirrorPOD) == 40, "layout mismatch");
  CUDA_CHECK(cudaMemcpy(&m, a.get_devObj(), sizeof(m), cudaMemcpyDeviceToHost));
  std::printf("mirror: _data=%p _start=%p _shift=%p\n", (void*)m._data,
              (void*)m._start, (void*)m._shift);

  std::ptrdiff_t hshift = 0;
  CUDA_CHECK(cudaMemcpy(&hshift, m._shift, sizeof(hshift), cudaMemcpyDeviceToHost));
  std::printf("device shift = %td\n", hshift);

  T* base = m._data;

  auto dump = [&](const char* tag, std::size_t begin, std::size_t cnt) {
    std::vector<T> h(cnt);
    CUDA_CHECK(cudaMemcpy(h.data(), base + begin, cnt * sizeof(T),
                          cudaMemcpyDeviceToHost));
    std::printf("  raw[%4zu..%4zu]:", begin, begin + cnt - 1);
    for (std::size_t k = 0; k < cnt; ++k) std::printf(" %.1f", (double)h[k]);
    std::printf("\n");
  };

  std::printf("== after upload (before rotates) ==\n");
  dump("A-head", 0, 6);
  dump("A-697", 695, 8);
  dump("A-tail", 995, 8);
  dump("B-head", 1023, 8);

  rotateLoopKernel<<<1, 1>>>(a.get_devObj(), 101);
  CUDA_CHECK(cudaDeviceSynchronize());

  std::printf("== after 101 rotates (offset 3) ==\n");
  CUDA_CHECK(cudaMemcpy(&hshift, m._shift, sizeof(hshift), cudaMemcpyDeviceToHost));
  std::printf("device shift = %td\n", hshift);
  dump("A-head", 0, 6);
  dump("A-697", 695, 8);
  dump("A-tail", 995, 8);
  dump("B-head", 1023, 8);

  // the view itself
  T* out = nullptr;
  CUDA_CHECK(cudaMalloc(&out, n * sizeof(T)));
  dumpViewKernel<<<2, 512>>>(a.get_devObj(), n, out);
  CUDA_CHECK(cudaDeviceSynchronize());
  std::vector<T> v(n);
  CUDA_CHECK(cudaMemcpy(v.data(), out, n * sizeof(T), cudaMemcpyDeviceToHost));
  std::printf("view[0..7]  :");
  for (int i = 0; i < 8; ++i) std::printf(" %.1f", (double)v[i]);
  std::printf("\nview[299..306]:");
  for (int i = 299; i < 307; ++i) std::printf(" %.1f", (double)v[i]);
  std::printf("\n");

  // reference: rotate by -303 over period 1024 -> view[i] = content[(721+i) mod 1024]
  std::printf("reference (period 1024) view[0..2]: 722 723 724\n");
  std::printf("reference (period 1000) view[0..2]: 698 699 1000->1\n");
  return 0;
}
