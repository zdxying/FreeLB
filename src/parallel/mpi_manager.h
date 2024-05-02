
/* This file is part of FreeLB, modified from OpenLB's mpiManager.h, with the following
 * copyright notice:
 *
 * // start of the original OpenLB's copyright notice
 *
 * This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 The OpenLB project
 *
 * // end of the original OpenLB's copyright notice
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

// mpi_manager.h

#pragma once

#ifdef MPI_ENABLED

#include <mpi.h>

#include <vector>
// std::unique_ptr
#include <memory>

template <typename T, unsigned int D>
class Vector;

#endif  // #ifdef MPI_ENABLED

#include <array>
#include <iostream>
#include <string>
// usleep
#include <unistd.h>

#ifdef MPI_ENABLED

#define MPI_RANK(x)           \
  if (mpi().getRank() != x) { \
    return;                   \
  }
// followed by a enclosing brace
#define IF_MPI_RANK(x) if (mpi().getRank() == x)

#ifdef MPI_DEBUG
#define MPI_DEBUG_WAIT        \
  if (mpi().getRank() == 0) { \
    int i = 1;                \
    while (i == 1) {          \
      sleep(1);               \
    }                         \
  }                           \
  mpi().barrier();
#else
#define MPI_DEBUG_WAIT
#endif  // #ifdef MPI_DEBUG

#else  // #ifdef MPI_ENABLED

#define MPI_RANK(x)
#define IF_MPI_RANK(x)
#define MPI_DEBUG_WAIT

#endif  // #ifdef MPI_ENABLED


#ifdef MPI_ENABLED

class MpiNonBlockingHelper {
 private:
  /// Size of the vector _mpiRequest/_mpiStatus
  unsigned _size;
  /// vector of MPI_Request
  std::unique_ptr<MPI_Request[]> _mpiRequest;
  /// vector of MPI_Status
  std::unique_ptr<MPI_Status[]> _mpiStatus;

 public:
  MpiNonBlockingHelper() : _size(0) {}
  ~MpiNonBlockingHelper() = default;

  MpiNonBlockingHelper(MpiNonBlockingHelper&& rhs) = default;
  MpiNonBlockingHelper(const MpiNonBlockingHelper&) = delete;
  MpiNonBlockingHelper& operator=(const MpiNonBlockingHelper&) = delete;

  /// Allocates memory
  void allocate(unsigned i);
  /// Reset
  void free() { _size = 0; }

  /// Returns the size of the vector _mpiRequest/_mpiStatus
  unsigned get_size() const { return _size; }

  /// Get the specified request object
  MPI_Request* get_mpiRequest(int i = 0) const;
  /// Get the specified status object
  MPI_Status* get_mpiStatus(int i = 0) const;

  void start(int i) { MPI_Start(get_mpiRequest(i)); }
  void wait(int i) { MPI_Wait(get_mpiRequest(i), get_mpiStatus(i)); }
  bool isDone(int i);

  /// Swap method
  void swap(MpiNonBlockingHelper& rhs);
};

class MpiManager {
 private:
  int numTasks, taskId;
  bool ok;

  friend MpiManager& mpi();

 public:
  MpiManager() : ok(false) {}
  ~MpiManager() {
    if (ok) {
      MPI_Finalize();
      ok = false;
    }
  }
  /// Initializes the mpi manager
  void init(int* argc, char*** argv);
  /// Returns the number of processes
  int getSize() const { return numTasks; }
  /// Returns the process ID
  int getRank() const { return taskId; }
  /// Returns process ID of main processor
  int bossId() const { return 0; }
  /// Tells whether current processor is main processor
  bool isMainProcessor() const { return bossId() == getRank(); }
  /// Returns universal MPI-time in seconds
  double getTime() const {
    if (!ok) return 0.;
    return MPI_Wtime();
  }

  /// Synchronizes the processes
  void barrier(MPI_Comm comm = MPI_COMM_WORLD) {
    if (!ok) return;
    MPI_Barrier(comm);
  }

  /// Synchronizes the processes and wait to ensure correct cout order
  void synchronizeIO(unsigned tDelay = 100, MPI_Comm comm = MPI_COMM_WORLD) {
    usleep(tDelay);
    barrier(comm);
  }

  /// Sends data at *buf, blocking
  template <typename T>
  void send(T* buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  template <typename... args>
  void send(std::vector<args...>& vec, int dest, int tag = 0,
            MPI_Comm comm = MPI_COMM_WORLD) {
    send(vec.data(), vec.size(), dest, tag, comm);
  }
  template <class T, std::size_t N>
  void send(std::array<T, N>& array, int dest, int tag = 0,
            MPI_Comm comm = MPI_COMM_WORLD) {
    send(array.data(), array.size(), dest, tag, comm);
  }

  /// Initialize persistent non-blocking send
  template <typename T>
  void sendInit(T* buf, int count, int dest, MPI_Request* request, int tag = 0,
                MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data at *buf, non blocking
  template <typename T>
  void iSend(T* buf, int count, int dest, MPI_Request* request, int tag = 0,
             MPI_Comm comm = MPI_COMM_WORLD);

  template <typename T, unsigned int D>
  void iSend(Vector<T, D>* buf, int arrsize, int dest, MPI_Request* request, int tag = 0,
             MPI_Comm comm = MPI_COMM_WORLD) {
    iSend(buf->data(), arrsize * D, dest, request, tag, comm);
  }

  /// Sends data at *buf, non blocking and buffered
  template <typename T>
  void ibSend(T* buf, int count, int dest, MPI_Request* request, int tag = 0,
              MPI_Comm comm = MPI_COMM_WORLD);

  /// Probe size of incoming message
  std::size_t probeReceiveSize(int source, MPI_Datatype type, int tag = 0,
                               MPI_Comm comm = MPI_COMM_WORLD);
  /// Probe size of incoming message with TYPE
  template <typename TYPE>
  std::size_t probeReceiveSize(int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Receives data at *buf, blocking
  template <typename T>
  void receive(T* buf, int count, int source, int tag = 0,
               MPI_Comm comm = MPI_COMM_WORLD);

  template <typename... args>
  void receive(std::vector<args...>& vec, int source, int tag = 0,
               MPI_Comm comm = MPI_COMM_WORLD) {
    receive(vec.data(), vec.size(), source, tag, comm);
  }
  template <class T, std::size_t N>
  void receive(std::array<T, N>& array, int source, int tag = 0,
               MPI_Comm comm = MPI_COMM_WORLD) {
    receive(array.data(), array.size(), source, tag, comm);
  }

  /// Initialize persistent non-blocking receive
  template <typename T>
  void recvInit(T* buf, int count, int dest, MPI_Request* request, int tag = 0,
                MPI_Comm comm = MPI_COMM_WORLD);

  /// Receives data at *buf, non blocking
  template <typename T>
  void iRecv(T* buf, int count, int source, MPI_Request* request, int tag = 0,
             MPI_Comm comm = MPI_COMM_WORLD);
  
  template <typename T, unsigned int D>
  void iRecv(Vector<T, D>* buf, int arrsize, int source, MPI_Request* request, int tag = 0,
             MPI_Comm comm = MPI_COMM_WORLD) {
    iRecv(buf->data(), arrsize * D, source, request, tag, comm);
  }

  /// Send and receive data between two partners
  template <typename T>
  void sendRecv(T* sendBuf, T* recvBuf, int count, int dest, int source, int tag = 0,
                MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data to master processor
  template <typename T>
  void sendToMaster(T* sendBuf, int sendCount, bool iAmRoot,
                    MPI_Comm comm = MPI_COMM_WORLD);

  /// Scatter data from one processor over multiple processors
  template <typename T>
  void scatterv(T* sendBuf, int* sendCounts, int* displs, T* recvBuf, int recvCount,
                int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Gather data from multiple processors to one processor
  template <typename T>
  void gatherv(T* sendBuf, int sendCount, T* recvBuf, int* recvCounts, int* displs,
               int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Broadcast data from one processor to multiple processors
  template <typename T>
  void bCast(T* sendBuf, int sendCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T>
  void bCast(T& sendVal, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Broadcast data when root is unknown to other processors
  template <typename T>
  void bCastThroughMaster(T* sendBuf, int sendCount, bool iAmRoot,
                          MPI_Comm comm = MPI_COMM_WORLD);

  /// Special case for broadcasting strings. Memory handling is automatic.
  void bCast(std::string& message, int root = 0);

  /// Reduction operation toward one processor
  template <typename T>
  void reduce(T& sendVal, T& recvVal, MPI_Op op, int root = 0, MPI_Comm = MPI_COMM_WORLD);

  /// Element-per-element reduction of a vector of data
  template <typename T>
  void reduceVect(std::vector<T>& sendVal, std::vector<T>& recvVal, MPI_Op op,
                  int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Reduction operation, followed by a broadcast
  template <typename T>
  void reduceAndBcast(T& reductVal, MPI_Op op, int root = 0,
                      MPI_Comm comm = MPI_COMM_WORLD);

  /// Complete a non-blocking MPI operation
  void wait(MPI_Request* request, MPI_Status* status) {
    if (!ok) return;
    MPI_Wait(request, status);
  }
  /// Complete a series of non-blocking MPI operations
  void waitAll(MpiNonBlockingHelper& mpiNbHelper) {
    if (!ok || mpiNbHelper.get_size() == 0) return;
    MPI_Waitall(mpiNbHelper.get_size(), mpiNbHelper.get_mpiRequest(),
                mpiNbHelper.get_mpiStatus());
  }
};


#else

class MpiManager {
 public:
  /// Initializes the mpi manager
  void init(int* argc, char*** argv) {}
  /// Returns the number of processes
  int getSize() const { return 1; }
  /// Returns the process ID
  int getRank() const { return 0; }
  /// Returns process ID of main processor
  int bossId() const { return 0; }
  /// Tells whether current processor is main processor
  bool isMainProcessor() const { return true; }

  /// Synchronizes the processes
  void barrier() const {};

  friend MpiManager& mpi();
};

#endif

MpiManager& mpi();