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

// field_statistics.h

#pragma once

#include "data_struct/field_struct.h"

template <typename BLOCKFIELDMANAGER>
class FieldStatistics {
private:
	const BLOCKFIELDMANAGER& BFM;

public:
	using FieldDataType = typename BLOCKFIELDMANAGER::datatype;
	using FloatType = typename BLOCKFIELDMANAGER::float_type;
	static constexpr unsigned int dim = BLOCKFIELDMANAGER::dim;

	FieldStatistics(const BLOCKFIELDMANAGER& bfm) : BFM(bfm) {}

	// get the average value of a field
	FieldDataType getAverage(std::size_t fieldidx = 0, bool useolap = true) const {
		FieldDataType TotalSum{};
		mpi().barrier();

		if (useolap) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalSum)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();

				FieldDataType sum{};
				std::size_t size{};
				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								sum += field[id];
								++id;
						}
					}
					size = (Nx - 2 * overlap) * (Ny - 2 * overlap);
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								sum += field[id];
								++id;
							}
						}
					}
					size = (Nx - 2 * overlap) * (Ny - 2 * overlap) * (Nz - 2 * overlap);
				}

				TotalSum += sum/size;
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalSum)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t size = field.size();
				FieldDataType sum{};
				for (std::size_t i = 0; i < size; ++i) {
					sum += field[i];
				}

				TotalSum += sum/size;
			}
		}

		TotalSum /= BFM.size();
#ifdef MPI_ENABLED
		mpi().reduceAndBcast(TotalSum, MPI_SUM);
		TotalSum /= mpi().getSize();
#endif
	return TotalSum;
	}

	// get the average value of a field
	template<typename FieldType, typename FlagType>
	FieldDataType getAverage(const BlockFieldManager<FieldType, FloatType, dim>& FlagFM, FlagType flag, std::size_t fieldidx = 0, bool useolap = true) const {
		FieldDataType TotalSum{};
		mpi().barrier();

		if (useolap) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalSum)
#endif
			for (std::size_t iField = 0; iField < BFM.getBlockFields().size(); ++iField) {
				const auto& blockfield = BFM.getBlockField(iField);
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();

				const auto& flagf = FlagFM.getBlockField(iField).getFieldType().getField(fieldidx);

				FieldDataType sum{};
				std::size_t size{};
				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (util::isFlag(flagf[id], flag)) {
									sum += field[id];
									++size;
								}
								++id;
						}
					}
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (util::isFlag(flagf[id], flag)) {
									sum += field[id];
									++size;
								}
								++id;
							}
						}
					}
				}

				FieldDataType temp = size == 0 ? FieldDataType{} : sum/size;
				TotalSum += temp;
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalSum)
#endif
			for (std::size_t iField = 0; iField < BFM.getBlockFields().size(); ++iField) {
				const auto& blockfield = BFM.getBlockField(iField);
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t fieldsize = field.size();

				const auto& flagf = FlagFM.getBlockField(iField).getFieldType().getField(fieldidx);

				std::size_t size{};
				FieldDataType sum{};
				for (std::size_t id = 0; id < fieldsize; ++id) {
					if (util::isFlag(flagf[id], flag)) {
						sum += field[id];
						++size;
					}
				}

				FieldDataType temp = size == 0 ? FieldDataType{} : sum/size;
				TotalSum += temp;
			}
		}

		TotalSum /= BFM.size();
#ifdef MPI_ENABLED
		mpi().reduceAndBcast(TotalSum, MPI_SUM);
		TotalSum /= mpi().getSize();
#endif
	return TotalSum;
	}

	// get the maximum value of a field
	FieldDataType getMax(std::size_t fieldidx = 0, bool useolap = true) const {
		FieldDataType TotalMax{};
		mpi().barrier();

		if (useolap) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : TotalMax)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();
				FieldDataType max{};

				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								max = max > field[id] ? max : field[id];
								++id;
						}
					}
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								max = max > field[id] ? max : field[id];
								++id;
							}
						}
					}
				}

				TotalMax = TotalMax > max ? TotalMax : max;
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : TotalMax)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t size = field.size();
				FieldDataType max{};
				for (std::size_t i = 0; i < size; ++i) {
					max = max > field[i] ? max : field[i];
				}

				TotalMax = TotalMax > max ? TotalMax : max;
			}
		}

#ifdef MPI_ENABLED
		mpi().reduceAndBcast(TotalMax, MPI_MAX);
#endif
		return TotalMax;
	}


	// get the minimum value of a field
	FieldDataType getMin(std::size_t fieldidx = 0, bool useolap = true) const {
		FieldDataType TotalMin{};
		mpi().barrier();

		if (useolap) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(min : TotalMin)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();
				FieldDataType min{};

				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								min = min < field[id] ? min : field[id];
								++id;
						}
					}
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								min = min < field[id] ? min : field[id];
								++id;
							}
						}
					}
				}

				TotalMin = TotalMin < min ? TotalMin : min;
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(min : TotalMin)
#endif			
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t size = field.size();
				FieldDataType min{};
				for (std::size_t i = 0; i < size; ++i) {
					min = min < field[i] ? min : field[i];
				}

				TotalMin = TotalMin < min ? TotalMin : min;
			}
		}

#ifdef MPI_ENABLED
		mpi().reduceAndBcast(TotalMin, MPI_MIN);
#endif
		return TotalMin;
	}


	// get number of a kind of value in a field
	std::size_t getCount(const FieldDataType value, std::size_t fieldidx = 0, bool useolap = true) const {
		std::size_t TotalValueCount{};
		mpi().barrier();

		if (useolap) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalValueCount)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();
				std::size_t count{};

				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (field[id] == value) ++count;
								++id;
						}
					}
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (field[id] == value) ++count;
								++id;
							}
						}
					}
				}

				TotalValueCount += count;
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalValueCount)
#endif			
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t size = field.size();
				std::size_t count{};
				for (std::size_t i = 0; i < size; ++i) {
					if (field[i] == value) ++count;
				}

				TotalValueCount += count;
			}
		}

#ifdef MPI_ENABLED
		mpi().reduceAndBcast(TotalValueCount, MPI_SUM);
#endif
		return TotalValueCount;
	}


	// get the percentage of a kind of value in a field
	FloatType getPercentage(const FieldDataType value, std::size_t fieldidx = 0, bool useolap = true) const {
		std::size_t TotalCount{};
		std::size_t TotalValueCount{};
		FloatType Result{};
		mpi().barrier();

		if (useolap) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalCount, TotalValueCount)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();
				std::size_t count{};

				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (field[id] == value) ++count;
								++id;
						}
					}
					TotalCount += (Nx - 2 * overlap) * (Ny - 2 * overlap);
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (field[id] == value) ++count;
								++id;
							}
						}
					}
					TotalCount += (Nx - 2 * overlap) * (Ny - 2 * overlap) * (Nz - 2 * overlap);
				}

				TotalValueCount += count;
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalCount, TotalValueCount)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t size = field.size();
				std::size_t count{};
				for (std::size_t i = 0; i < size; ++i) {
					if (field[i] == value) ++count;
				}

				TotalCount += size;
				TotalValueCount += count;
			}
		}

		Result = static_cast<FloatType>(TotalValueCount) / TotalCount;

#ifdef MPI_ENABLED
		mpi().reduceAndBcast(Result, MPI_SUM);
		Result /= mpi().getSize();
#endif
		return Result;
	}


	// get number of a kind of flag in a flag field
	std::size_t getFlagCount(const FieldDataType flag, std::size_t fieldidx = 0, bool useolap = true) const {
		std::size_t TotalValueCount{};
		mpi().barrier();

		if (useolap) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalValueCount)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();
				std::size_t count{};

				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (util::isFlag(field[id], flag)) ++count;
								++id;
						}
					}
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (util::isFlag(field[id], flag)) ++count;
								++id;
							}
						}
					}
				}

				TotalValueCount += count;
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalValueCount)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t size = field.size();
				std::size_t count{};
				for (std::size_t i = 0; i < size; ++i) {
					if (util::isFlag(field[i], flag)) ++count;
				}

				TotalValueCount += count;
			}
		}

#ifdef MPI_ENABLED
		mpi().reduceAndBcast(TotalValueCount, MPI_SUM);
#endif
		return TotalValueCount;
	}


	// get the percentage of a kind of flag value in a flag field
	FloatType getFlagPercentage(const FieldDataType flag, std::size_t fieldidx = 0, bool useolap = true) const {
		std::size_t TotalCount{};
		std::size_t TotalValueCount{};
		FloatType Result{};
		mpi().barrier();

		if (useolap) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalCount, TotalValueCount)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();
				std::size_t count{};

				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (util::isFlag(field[id], flag)) ++count;
								++id;
						}
					}
					TotalCount += (Nx - 2 * overlap) * (Ny - 2 * overlap);
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								if (util::isFlag(field[id], flag)) ++count;
								++id;
							}
						}
					}
					TotalCount += (Nx - 2 * overlap) * (Ny - 2 * overlap) * (Nz - 2 * overlap);
				}

				TotalValueCount += count;
			}
		} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalCount, TotalValueCount)
#endif
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t size = field.size();
				std::size_t count{};
				for (std::size_t i = 0; i < size; ++i) {
					if (util::isFlag(field[i], flag)) ++count;
				}

				TotalCount += size;
				TotalValueCount += count;
			}
		}

		Result = static_cast<FloatType>(TotalValueCount) / TotalCount;

#ifdef MPI_ENABLED
		mpi().reduceAndBcast(Result, MPI_SUM);
		Result /= mpi().getSize();
#endif
		return Result;
	}


	// print statistics of a flag field
	void printFlagStatistics(std::size_t fieldidx = 0, bool printpercentage = true, bool useolap = true) const {
		std::size_t TotalCount{};
		std::size_t Total{};
		// get the number of different flags
		std::vector<std::pair<FieldDataType, std::size_t>> flagcount;
		std::vector<std::pair<FieldDataType, std::size_t>> Result;
		mpi().barrier();

		if (useolap) {
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();

				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
					for (int j = overlap; j < Ny - overlap; ++j) {
						id = j * Nx + overlap;
						for (int i = overlap; i < Nx - overlap; ++i) {
							// find if the flag is in the vector
							const auto var = field[id];
							auto it = std::find_if(flagcount.begin(), flagcount.end(), [var](const auto& p) {
								return util::isFlag(var, p.first);
							});
							if (it != flagcount.end()) {
								++it->second;
							} else {
								flagcount.push_back({var, std::size_t{1}});
							}
							++id;
						}
					}
					TotalCount += (Nx - 2 * overlap) * (Ny - 2 * overlap);
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								// find if the flag is in the vector
								const auto var = field[id];
								auto it = std::find_if(flagcount.begin(), flagcount.end(), [var](const auto& p) {
									return util::isFlag(var, p.first);
								});
								if (it != flagcount.end()) {
									++it->second;
								} else {
									flagcount.push_back({var, std::size_t{1}});
								}
								++id;
							}
						}
					}
					TotalCount += (Nx - 2 * overlap) * (Ny - 2 * overlap) * (Nz - 2 * overlap);
				}
			}
		} else {
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t size = field.size();
				for (std::size_t i = 0; i < size; ++i) {
					// find if the flag is in the vector
					const auto& var = field[i];
					auto it = std::find_if(flagcount.begin(), flagcount.end(), [var](const auto& p) {
						return util::isFlag(var, p.first);
					});
					if (it != flagcount.end()) {
						++it->second;
					} else {
						flagcount.push_back({var, std::size_t{1}});
					}
				}
				TotalCount += size;
			}
		}

#ifdef MPI_ENABLED
		// merge all ranks' vector
		mpi().reduce(TotalCount, Total, MPI_SUM);
		// 1. prepare send buffers
		std::vector<FieldDataType> flagsendbuf;
		std::vector<std::size_t> countsendbuf;
		for (const auto& p : flagcount) {
			flagsendbuf.push_back(p.first);
			countsendbuf.push_back(p.second);
		} 
		// 2. determine the size of send buffers
		int sendbuf_size = static_cast<int>(flagsendbuf.size());
		// 3. gather the sizes at the main rank
		std::vector<int> recvcounts(mpi().getSize());
		MPI_Gather(&sendbuf_size, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		// 4. allocate space for receiving all sendbuf at the main rank
		std::vector<int> displs(mpi().getSize(), 0);
		int total_size{};
    IF_MPI_RANK(0) {
			for (int i = 0; i < mpi().getSize(); ++i) {
				displs[i] = total_size;
				total_size += recvcounts[i];
			}
    }
		std::vector<FieldDataType> flagrecvbuf(total_size);
		std::vector<std::size_t> countrecvbuf(total_size);
		// 5. gather all sendbuf at the main rank
		mpi().gatherv(flagsendbuf.data(), sendbuf_size, flagrecvbuf.data(), recvcounts.data(), displs.data());
		mpi().gatherv(countsendbuf.data(), sendbuf_size, countrecvbuf.data(), recvcounts.data(), displs.data());
		// 6. process at the main rank
		MPI_RANK(0);
		for (int i = 0; i < total_size; ++i) {
			const auto& var = flagrecvbuf[i];
			auto it = std::find_if(Result.begin(), Result.end(), [var](const auto& p) {
				return util::isFlag(var, p.first);
			});
			if (it != Result.end()) {
				it->second += countrecvbuf[i];
			} else {
				Result.push_back({var, countrecvbuf[i]});
			}
		}
#else
		Result = flagcount;
		Total = TotalCount;
#endif
		// sort the vector
		std::sort(Result.begin(), Result.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
		// print info
		std::cout << "[Flag statistics]:" << std::endl;
		std::cout << "  Flag | ";
		if (printpercentage) std::cout << "Percent% | ";
		std::cout << "Count" << std::endl;
		for (const auto& p : Result) {
			std::cout << "  " << std::setw(7) << std::left << static_cast<int>(p.first);
			if (printpercentage) std::cout << std::setw(11) << std::left << static_cast<FloatType>(p.second) / Total * 100;
			std::cout << p.second << std::endl;
		}
	}


	// print statistics of a field
	void printValueStatistics(std::size_t fieldidx = 0, bool printpercentage = true, bool useolap = true) const {
		std::size_t TotalCount{};
		std::size_t Total{};
		// get the number of different flags
		std::vector<std::pair<FieldDataType, std::size_t>> valuecount;
		std::vector<std::pair<FieldDataType, std::size_t>> Result;
		mpi().barrier();

		if (useolap) {
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const auto& blockxd = blockfield.getBlock();
				const int overlap = blockxd.getOverlap();
				const int Nx = blockxd.getNx();
				const int Ny = blockxd.getNy();
				const int Nz = blockxd.getNz();

				std::size_t id{};
				if constexpr (BLOCKFIELDMANAGER::dim == 2) {
					for (int j = overlap; j < Ny - overlap; ++j) {
						id = j * Nx + overlap;
						for (int i = overlap; i < Nx - overlap; ++i) {
							// find if the flag is in the vector
							const auto var = field[id];
							auto it = std::find_if(valuecount.begin(), valuecount.end(), [var](const auto& p) {
								return var == p.first;
							});
							if (it != valuecount.end()) {
								++it->second;
							} else {
								valuecount.push_back({var, std::size_t{1}});
							}
							++id;
						}
					}
					TotalCount += (Nx - 2 * overlap) * (Ny - 2 * overlap);
				} else if constexpr (BLOCKFIELDMANAGER::dim == 3) {
					std::size_t NxNy = Nx * Ny;
					for (int k = overlap; k < Nz - overlap; ++k) {
						for (int j = overlap; j < Ny - overlap; ++j) {
							id = k * NxNy + j * Nx + overlap;
							for (int i = overlap; i < Nx - overlap; ++i) {
								// find if the flag is in the vector
								const auto var = field[id];
								auto it = std::find_if(valuecount.begin(), valuecount.end(), [var](const auto& p) {
									return var == p.first;
								});
								if (it != valuecount.end()) {
									++it->second;
								} else {
									valuecount.push_back({var, std::size_t{1}});
								}
								++id;
							}
						}
					}
					TotalCount += (Nx - 2 * overlap) * (Ny - 2 * overlap) * (Nz - 2 * overlap);
				}
			}
		} else {
			for (const auto& blockfield : BFM.getBlockFields()) {
				const auto& field = blockfield.getFieldType().getField(fieldidx);
				const std::size_t size = field.size();
				for (std::size_t i = 0; i < size; ++i) {
					// find if the flag is in the vector
					const auto& var = field[i];
					auto it = std::find_if(valuecount.begin(), valuecount.end(), [var](const auto& p) {
						return var == p.first;
					});
					if (it != valuecount.end()) {
						++it->second;
					} else {
						valuecount.push_back({var, std::size_t{1}});
					}
				}
				TotalCount += size;
			}
		}
		
#ifdef MPI_ENABLED
		// merge all ranks' vector
		mpi().reduce(TotalCount, Total, MPI_SUM);
		// 1. prepare send buffers
		std::vector<FieldDataType> valuesendbuf;
		std::vector<std::size_t> countsendbuf;
		for (const auto& p : valuecount) {
			valuesendbuf.push_back(p.first);
			countsendbuf.push_back(p.second);
		} 
		// 2. determine the size of send buffers
		int sendbuf_size = static_cast<int>(valuesendbuf.size());
		// 3. gather the sizes at the main rank
		std::vector<int> recvcounts(mpi().getSize());
		MPI_Gather(&sendbuf_size, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		// 4. allocate space for receiving all sendbuf at the main rank
		std::vector<int> displs(mpi().getSize(), 0);
		int total_size{};
    IF_MPI_RANK(0) {
			for (int i = 0; i < mpi().getSize(); ++i) {
				displs[i] = total_size;
				total_size += recvcounts[i];
			}
    }
		std::vector<FieldDataType> valuerecvbuf(total_size);
		std::vector<std::size_t> countrecvbuf(total_size);
		// 5. gather all sendbuf at the main rank
		mpi().gatherv(valuesendbuf.data(), sendbuf_size, valuerecvbuf.data(), recvcounts.data(), displs.data());
		mpi().gatherv(countsendbuf.data(), sendbuf_size, countrecvbuf.data(), recvcounts.data(), displs.data());
		// 6. process at the main rank
		MPI_RANK(0);
		for (int i = 0; i < total_size; ++i) {
			const auto& var = valuerecvbuf[i];
			auto it = std::find_if(Result.begin(), Result.end(), [var](const auto& p) {
				return var == p.first;
			});
			if (it != Result.end()) {
				it->second += countrecvbuf[i];
			} else {
				Result.push_back({var, countrecvbuf[i]});
			}
		}
#else
		Result = valuecount;
		Total = TotalCount;
#endif
		// sort the vector
		std::sort(Result.begin(), Result.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
		// print info
		std::cout << "[Value statistics]:" << std::endl;
		std::cout << "Value | ";
		if (printpercentage) std::cout << "Percent% | ";
		std::cout << "Count" << std::endl;
		for (const auto& p : Result) {
			std::cout << std::setw(8) << std::left << p.first;
			if (printpercentage) std::cout << std::setw(11) << std::left << static_cast<FloatType>(p.second) / Total * 100;
			std::cout << p.second << std::endl;
		}
	}

};