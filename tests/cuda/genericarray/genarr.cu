

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;

std::size_t N = 10;


// GenericArray
template <typename T>
__any__ void addGenericArrayImp(cudev::GenericArray<T> &a, std::size_t id, T value) {
  a[id] += value;
}
template <typename T>
__global__ void addGenericArray_kernel(cudev::GenericArray<T> *a, T value) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < a->size()) {
    addGenericArrayImp(*a, idx, value);
  }
}
template <typename T>
void addGenericArray(GenericArray<T> &a, T value) {
  const unsigned int blockSize = 32;
  const unsigned int blockNum = (a.size() + blockSize - 1) / blockSize;
  addGenericArray_kernel<<<blockNum, blockSize>>>(a.get_devObj(), value);
}

// Data
template <typename T, typename Base>
void addData(Data<T, Base> &a, T value) {
  addData_kernel<<<1, 1>>>(a.get_devObj(), value);
}
template <typename T, typename Base>
__global__ void addData_kernel(cudev::Data<T, Base> *a, T value) {
  addDataImp(*a, value);
}
template <typename T, typename Base>
__any__ void addDataImp(cudev::Data<T, Base> &a, T value) {
  a.get() += value;
}

// Genericvector
template <typename T>
__any__ void addvectorImp(cudev::Genericvector<T> &a, std::size_t id, T value) {
  a[id] += value;
}
template <typename T>
__global__ void addvector_kernel(cudev::Genericvector<T> *a, T value) {
  std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < a->size()) {
    addvectorImp(*a, idx, value);
  }
}
template <typename T>
void addvector(Genericvector<T> &a, T value) {
  const unsigned int blockSize = 32;
  const unsigned int blockNum = (a.size() + blockSize - 1) / blockSize;
  addvector_kernel<<<blockNum, blockSize>>>(a.get_devObj(), value);
}

// StreamArray
template <typename T>
__any__ void addStreamArrayImp(cudev::StreamArray<T> &a, std::size_t id) {
  a[id] = id + 1;
}
template <typename T>
__global__ void addStreamArray_kernel(cudev::StreamArray<T> *a) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < a->size()) {
    addStreamArrayImp(*a, idx);
  }
}
template <typename T>
void addStreamArray(StreamArray<T> &a) {
  const unsigned int blockSize = 32;
  const unsigned int blockNum = (a.size() + blockSize - 1) / blockSize;
  addStreamArray_kernel<<<blockNum, blockSize>>>(a.get_devObj());
}

void set(StreamArray<T> &arr) {
  for(int i = 0; i < arr.size(); i++) {
    arr[i] = i;
  } 
}
// void print(StreamArray<T> &arr) {
//   for(int i = 0; i < arr.size(); i++) {
//     std::cout << arr[i] << " ";
//   }
//   std::cout << std::endl;
// }

// GenericField
template <typename ArrayType, typename T, typename Base>
__any__ void addGenericFieldImp(cudev::GenericField<ArrayType, Base> &a, std::size_t id, T value) {
  a.get(id) += value;
}
template <typename ArrayType, typename T, typename Base>
__global__ void addGenericField_kernel(cudev::GenericField<ArrayType, Base> *a, T value) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < a->getField().size()) {
    addGenericFieldImp(*a, idx, value);
  }
}
template <typename ArrayType, typename T, typename Base>
void addGenericField(GenericField<ArrayType, Base> &a, T value) {
  const unsigned int blockSize = 32;
  const unsigned int blockNum = (a.getField().size() + blockSize - 1) / blockSize;
  addGenericField_kernel<<<blockNum, blockSize>>>(a.get_devObj(), value);
}

int main() {
  // std::cout << "sizeof(T) = " << sizeof(T) << std::endl;
  // std::cout << "sizeof(Vector<T, 3>) = " << sizeof(Vector<T, 3>) << std::endl;

  VELOCITY<T, 3> v(N, Vector<T, 3>{T{1.}, T{2.}, T{3.}});
  std::cout << v.get(N-1)[0] << " " << v.get(N-1)[1] << " " << v.get(N-1)[2] << std::endl;
  addGenericField(v, Vector<T, 3>{T{1.}, T{2.}, T{3.}});
  v.copyToHost();
  std::cout << v.get(N-1)[0] << " " << v.get(N-1)[1] << " " << v.get(N-1)[2] << std::endl;

  // RHO<T> v(N, T{1.});
  // CONSTRHO<T> v(1, T{1.});
  // addData(v, T{1.});
  // v.copyToHost();
  // std::cout << v.get() << std::endl;

  // Genericvector<int> v(N, 1);
  // std::cout << v[N-1] << std::endl;
  // addvector(v, 1);
  // v.copyToHost();
  // std::cout << v[N-1] << std::endl;

  // StreamArray<T> sarr(100000, T{});
  // int offset = int(std::sqrt(100000));
  // sarr.setOffset(offset);
  // set(sarr);
  // sarr.copyToDevice();
  // // print(sarr);
  // addStreamArray(sarr, T{1.});
  // sarr.dev_rotate();
  // sarr.copyToHost();
  // print(sarr);  
  // Timer MainLoopTimer;
  // MainLoopTimer.START_TIMER();
  // for(int i = 0; i < 10000; ++i) {
  //   addStreamArray(sarr);
  //   // sarr.dev_rotate();
  //   sarr.rotate_dev();
  // }
  // MainLoopTimer.END_TIMER();
  // std::cout << MainLoopTimer.GetDurationCount_Only() <<  std::endl;


  return 0;
}

// ----------------------------------------------------------------

// __global__ void childKernel() {
//     printf("Hello from child kernel\n");
// }

// __global__ void parentKernel() {
//     // Launch child kernel from within the parent kernel
//     childKernel<<<2, 1>>>();
//     printf("Hello from parent kernel\n");
// }

// int main() {
//     // Launch the parent kernel
//     parentKernel<<<2, 1>>>();
//     // Wait for the parent kernel to finish
//     cudaDeviceSynchronize();
//     return 0;
// }

// ----------------------------------------------------------------

// class vecarray {
//  public:
//   int *vecptr;  // array of pointers pointing to array
//   int dim;      // store length of each array pointed to

//   __device__ __host__ vecarray();  // constructor
//   __device__ __host__ int sum();   // sum up all the elements in the array being
//                                   // pointed to
// };

// vecarray::vecarray() {
//   vecptr = NULL;
//   dim = 0;
// }

// __device__ __host__ int vecarray::sum() {
//   int j = 0, s = 0;
//   for (j = 0; j < dim; j++) s += vecptr[j];
//   return s;
// }

// __global__ void addvecarray(vecarray *v, int *s) { *s = v->sum(); }

// int main() {        // copy *V to device, do sum() and pass back
//   vecarray *dev_v;  // the result by dev_v
//   vecarray v;
//   int a[3] = {1, 2, 3};  // initialize v manually
//   int result = 0;
//   int *dev_result;
//   v.vecptr = a;
//   v.dim = 3;
//   int *vptr;

//   cudaMalloc((void **)&dev_v, sizeof(vecarray));
//   cudaMemcpy(dev_v, &v, sizeof(vecarray), cudaMemcpyHostToDevice);  // copy class
//   object

// 	cudaMalloc((void **)&(vptr), v.dim * sizeof(int));
//   // copy arrays
//   cudaMemcpy(vptr, v.vecptr, v.dim * sizeof(int), cudaMemcpyHostToDevice);

//   cudaMemcpy(&(dev_v->vecptr), &vptr, sizeof(int *), cudaMemcpyHostToDevice);


//   cudaMalloc((void **)&dev_result, sizeof(int));
//   addvecarray<<<1, 1>>>(dev_v, dev_result);

//   cudaMemcpy(&result, dev_result, sizeof(int), cudaMemcpyDeviceToHost);
//   printf("the result is %d\n", result);
//   return 0;
// }