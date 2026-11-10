// Probe v2: cost of splitting one kernel launch into N launches (GPU multi-block proxy)
// Adds warmup, CUDA-Graph comparison, and variance reporting.

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>

#define CHECK(e) do { cudaError_t _e = (e); if (_e != cudaSuccess) { \
  std::fprintf(stderr, "cuda error %s @ %d\n", cudaGetErrorString(_e), __LINE__); std::exit(1); } } while (0)

__global__ void copy_kernel(float* out, const float* in, std::size_t n) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) out[i] = in[i] * 2.0f;
}

__global__ void empty_kernel() {}

static cudaStream_t S = nullptr;

static double run(int N, int reps, bool use_graph, bool real_work,
                  float* dout, const float* din, std::size_t total) {
  const unsigned int bs = 128;
  const std::size_t per = (total + N - 1) / N;
  const unsigned int grid = (unsigned int)((per + bs - 1) / bs);

  cudaGraph_t graph = nullptr;
  cudaGraphExec_t exec = nullptr;
  if (use_graph) {
    CHECK(cudaStreamSynchronize(S));
    CHECK(cudaStreamBeginCapture(S, cudaStreamCaptureModeRelaxed));
    for (int k = 0; k < N; ++k) {
      if (real_work) {
        std::size_t off = (std::size_t)k * per;
        std::size_t n = (off + per > total) ? (total - off) : per;
        copy_kernel<<<grid, bs, 0, S>>>(dout + off, din + off, n);
      } else {
        empty_kernel<<<1, 32, 0, S>>>();
      }
    }
    CHECK(cudaStreamEndCapture(S, &graph));
    CHECK(cudaGraphInstantiate(&exec, graph, nullptr, nullptr, 0));
  }

  auto body = [&]() {
    if (use_graph) {
      CHECK(cudaGraphLaunch(exec, S));
    } else {
      for (int k = 0; k < N; ++k) {
        if (real_work) {
          std::size_t off = (std::size_t)k * per;
          std::size_t n = (off + per > total) ? (total - off) : per;
          copy_kernel<<<grid, bs, 0, S>>>(dout + off, din + off, n);
        } else {
          empty_kernel<<<1, 32, 0, S>>>();
        }
      }
    }
  };

  for (int r = 0; r < 20; ++r) body();  // warmup
  CHECK(cudaDeviceSynchronize());

  cudaEvent_t a, b;
  CHECK(cudaEventCreate(&a));
  CHECK(cudaEventCreate(&b));
  CHECK(cudaEventRecord(a));
  for (int r = 0; r < reps; ++r) body();
  CHECK(cudaEventRecord(b));
  CHECK(cudaEventSynchronize(b));
  float ms = 0.f;
  CHECK(cudaEventElapsedTime(&ms, a, b));
  CHECK(cudaEventDestroy(a));
  CHECK(cudaEventDestroy(b));

  if (use_graph) {
    CHECK(cudaGraphExecDestroy(exec));
    CHECK(cudaGraphDestroy(graph));
  }
  return (double)ms * 1000.0 / reps;  // us per step
}

int main() {
  int dev = 0;
  CHECK(cudaGetDevice(&dev));
  CHECK(cudaStreamCreate(&S));
  int sm = 0, regs = 0, maxthr = 0;
  cudaDeviceGetAttribute(&sm, cudaDevAttrMultiProcessorCount, dev);
  cudaDeviceGetAttribute(&regs, cudaDevAttrMaxRegistersPerMultiprocessor, dev);
  cudaDeviceGetAttribute(&maxthr, cudaDevAttrMaxThreadsPerMultiProcessor, dev);
  cudaDeviceProp p{};
  CHECK(cudaGetDeviceProperties(&p, dev));
  std::printf("device=%s  SM=%d  regs/SM=%d  maxThreads/SM=%d  clock=%.2f GHz\n",
              p.name, sm, regs, maxthr, p.clockRate / 1e6);

  const std::size_t total = 1u << 22;  // 4M floats
  float *din, *dout;
  CHECK(cudaMalloc(&din, total * sizeof(float)));
  CHECK(cudaMalloc(&dout, total * sizeof(float)));
  CHECK(cudaMemset(din, 1, total * sizeof(float)));

  int Ns[] = {1, 2, 4, 8, 16, 32, 64, 128};

  std::printf("\n[A] empty kernel, us/step (pure launch cost)\n");
  std::printf("%6s %12s %12s %12s %10s\n", "N", "plain", "graph", "ns/launch", "graph eff");
  for (int N : Ns) {
    double tp = run(N, 500, false, false, dout, din, total);
    double tg = run(N, 500, true, false, dout, din, total);
    std::printf("%6d %12.2f %12.2f %12.1f %9.2fx\n", N, tp, tg, tp * 1000.0 / N, tp / tg);
  }

  std::printf("\n[B] DRAM-bound work split into N kernels, total=%zu elems\n", total);
  std::printf("%6s %10s %10s %8s %10s %10s %8s\n",
              "N", "us/step", "GB/s", "vs N=1", "graph us", "graph GB/s", "gain");
  double base = 0.0;
  for (int N : Ns) {
    double tp = run(N, 100, false, true, dout, din, total);
    double tg = run(N, 100, true, true, dout, din, total);
    if (N == 1) base = tp;
    double gbs = (double)total * 2.0 * sizeof(float) / (tp * 1e-6) / 1e9;
    double gbsg = (double)total * 2.0 * sizeof(float) / (tg * 1e-6) / 1e9;
    std::printf("%6d %10.2f %10.1f %7.2fx %10.2f %10.1f %7.2fx\n",
                N, tp, gbs, tp / base, tg, gbsg, tp / tg);
  }

  CHECK(cudaFree(din));
  CHECK(cudaFree(dout));
  return 0;
}
