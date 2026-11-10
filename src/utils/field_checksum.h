// field_checksum.h -- pop-field checksum (sum / sumsq / maxabs) for solver
// regression and cross-implementation comparison.
//
// sum   = sum of all distribution values       ~ total fluid mass (N * rho)
// sumsq = sum of squared values                sensitive fingerprint of the
//                                              field's fine structure
// maxabs = largest |f|                          boundedness / blow-up sentinel
//
// The GPU and CPU variants accumulate in the same per-cell order with double
// accumulators, so two solvers that agree bitwise produce identical checksums
// and the difference between two runs quantifies field divergence.
//
// Include after freelb.h (needs the POP alias and the block-lattice types).

#pragma once

#include <array>
#include <cstddef>
#include <iomanip>
#include <iostream>

namespace frelb_diag {

#ifdef __CUDACC__

template <typename T, typename LatSet, typename TypePack>
__global__ void PopChecksumKernel(
    cudev::BlockLattice<T, LatSet, TypePack>* blocklat, std::size_t n,
    double* out) {
  double s1 = 0.0, s2 = 0.0, mx = 0.0, nanc = 0.0;
  long long firstnan = -1;
  // inside the device, POP must resolve to the cudev alias
  auto& popf = blocklat->template getField<cudev::POP<T, LatSet::q>>();
  for (std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += (std::size_t)gridDim.x * blockDim.x) {
    for (unsigned int d = 0; d < LatSet::q; ++d) {
      const double v = static_cast<double>(popf.getField(d).getdataPtr(i)[0]);
      if (v != v) {  // NaN sentinel: count + first (cell, pop) occurrence
        ++nanc;
        if (firstnan < 0) firstnan = (long long)(i * (std::size_t)LatSet::q + d);
        continue;
      }
      s1 += v;
      s2 += v * v;
      const double a = v < 0 ? -v : v;
      mx = a > mx ? a : mx;
    }
  }
  atomicAdd(&out[0], s1);
  atomicAdd(&out[1], s2);
  // non-negative doubles, so an integer atomicMax on the bit pattern works
  atomicMax(reinterpret_cast<unsigned long long*>(&out[2]),
            __double_as_longlong(mx));
  atomicAdd(&out[3], nanc);
  atomicMin(reinterpret_cast<unsigned long long*>(&out[4]),
            firstnan < 0 ? 0x7fffffffffffffffLL : (unsigned long long)firstnan);
}

// device pop-field checksum; the POP field must already be on the device
// NaN diagnostics from the last PopChecksumDevice call
inline double nan_count_ = 0.0;
inline long long first_nan_ = -1;

template <typename BLOCKLAT>
std::array<double, 3> PopChecksumDevice(BLOCKLAT& bl) {
  using T = typename BLOCKLAT::FloatType;
  using LS = typename BLOCKLAT::LatticeSet;
  using TP = typename BLOCKLAT::cudev_TypePack;
  double* d_stats = nullptr;
  cudaMalloc(&d_stats, 5 * sizeof(double));
  cudaMemset(d_stats, 0, 5 * sizeof(double));
  const std::size_t n = bl.getN();
  PopChecksumKernel<T, LS, TP><<<1024, 256>>>(bl.get_devObj(), n, d_stats);
  std::array<double, 5> h5{0, 0, 0, 0, 0};
  cudaMemcpy(h5.data(), d_stats, 5 * sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(d_stats);
  nan_count_ = h5[3];
  first_nan_ = (long long)h5[4];
  return {h5[0], h5[1], h5[2]};
}

#endif  // __CUDACC__

// host pop-field checksum; popf is the block lattice's POP GenericField,
// q the lattice direction count
template <typename POPFIELD>
std::array<double, 3> PopChecksumHost(POPFIELD& popf, std::size_t n,
                                      unsigned int q) {
  std::array<double, 3> s{0, 0, 0};
  for (std::size_t i = 0; i < n; ++i) {
    for (unsigned int d = 0; d < q; ++d) {
      const double v = static_cast<double>(popf.getField(d)[i]);
      s[0] += v;
      s[1] += v * v;
      const double a = v < 0 ? -v : v;
      if (a > s[2]) s[2] = a;
    }
  }
  return s;
}

inline void PrintChecksum(const char* tag, const std::array<double, 3>& s) {
  std::cout << tag << " checksum: sum " << std::setprecision(10) << s[0]
            << "  sumsq " << s[1] << "  maxabs " << s[2] << std::endl;
}

}  // namespace frelb_diag
