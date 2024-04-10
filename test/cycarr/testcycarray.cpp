#include "data_struct/field.h"

template <typename T>
void printcyc(CyclicArray<T> &arr, int size) {
  for (int i = 0; i < size; ++i) {
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;
}

int main() {
  int size = 10;
  CyclicArray<int> arr(size, 0);
  for (int i = 0; i < size; ++i) {
    arr[i] = arr[i] + i;
  }
  printcyc(arr, size);
  // std::cout << arr[-5]<<std::endl;
  arr.rotate(long(-3));
  printcyc(arr, size);
	int x = arr.getPrevious(6);
  std::cout << x << std::endl;
  return 0;
}