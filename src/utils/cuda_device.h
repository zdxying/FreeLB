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

#pragma once

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <cstdint>
#include <iostream>
#include <string>

// number of threads per block
#define THREADS_PER_BLOCK 32

template <typename T, unsigned int D>
class Vector;

// check cuda status
void check_cuda_status(std::string message = "FreeLB cuda_device.h") {
  cudaError_t status = cudaGetLastError();
  if (status != cudaSuccess) {
    std::cerr << "[" << message << "] CUDA error: " << cudaGetErrorString(status)
              << std::endl;
    exit(1);
  }
}

// malloc on device
template <typename T>
T* cuda_malloc(std::size_t N) {
  T* ptr;
  cudaMalloc(&ptr, N * sizeof(T));
  check_cuda_status("cuda_malloc");
  return ptr;
}

// free on device
template <typename T>
void cuda_free(T* ptr) {
  cudaFree(ptr);
  check_cuda_status("cuda_free");
}

template <typename T>
__global__ void fillKernel(T* deviceData, T value, std::size_t size) {
  const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < size) {
    deviceData[idx] = value;
  }
}

template <typename T>
void cuda_fill(T* deviceData, std::size_t size, T value) {
  const unsigned int blockSize = THREADS_PER_BLOCK;
  const unsigned int blockNum = (size + blockSize - 1) / blockSize;
  fillKernel<<<blockNum, blockSize>>>(deviceData, value, size);
  cudaDeviceSynchronize();
  check_cuda_status("cuda_fill");
}

template <typename T>
void host_to_device(T* deviceData, T* hostData, std::size_t size = 1) {
  cudaMemcpy(deviceData, hostData, size * sizeof(T), cudaMemcpyHostToDevice);
  check_cuda_status("host_to_device");
}

template <typename T, unsigned int D>
void host_to_device(Vector<T, D>* deviceData, Vector<T, D>* hostData, std::size_t size = 1) {
  cudaMemcpy(deviceData->data(), hostData->data(), D * size * sizeof(T),
			 cudaMemcpyHostToDevice);
  check_cuda_status("Vector host_to_device");
}

template <typename T>
void device_to_host(T* hostData, T* deviceData, std::size_t size = 1) {
  cudaMemcpy(hostData, deviceData, size * sizeof(T), cudaMemcpyDeviceToHost);
  check_cuda_status("device_to_host");
}

template <typename T, unsigned int D>
void device_to_host(Vector<T, D>* hostData, Vector<T, D>* deviceData, std::size_t size = 1) {
  cudaMemcpy(hostData->data(), deviceData->data(), D * size * sizeof(T),
			 cudaMemcpyDeviceToHost);
  check_cuda_status("Vector device_to_host");
}

template <typename T>
void device_to_device(T* dst, T* src, std::size_t size = 1) {
  cudaMemcpy(dst, src, size * sizeof(T), cudaMemcpyDeviceToDevice);
  check_cuda_status("device_to_device");
}

template <typename T, unsigned int D>
void device_to_device(Vector<T, D>* dst, Vector<T, D>* src, std::size_t size = 1) {
  cudaMemcpy(dst->data(), src->data(), D * size * sizeof(T), cudaMemcpyDeviceToDevice);
  check_cuda_status("Vector device_to_device");
}

template <typename T>
__device__ void dev_copy(T* dst, T* src, std::size_t size) {
  for (std::size_t i = 0; i < size; ++i) {
    dst[i] = src[i];
  }
}

template <typename T>
__global__ void copyKernel(T* dst, T* src, std::size_t size) {
  const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < size) {
    dst[idx] = src[idx];
  }
}