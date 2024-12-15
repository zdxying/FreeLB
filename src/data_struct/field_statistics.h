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

	FieldStatistics(const BLOCKFIELDMANAGER& bfm) : BFM(bfm) {}

	// get the average value of a field
	FieldDataType getAverage(std::size_t fieldidx = 0, bool useolap = true) const {
		FieldDataType TotalSum{};

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

		return TotalSum / BFM.size();
	}


	// get the maximum value of a field
	FieldDataType getMax(std::size_t fieldidx = 0, bool useolap = true) const {
		FieldDataType TotalMax{};

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

		return TotalMax;
	}


	// get the minimum value of a field
	FieldDataType getMin(std::size_t fieldidx = 0, bool useolap = true) const {
		FieldDataType TotalMin{};

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

		return TotalMin;
	}


	// get number of a kind of value in a field
	std::size_t getCount(const FieldDataType value, std::size_t fieldidx = 0, bool useolap = true) const {
		std::size_t TotalValueCount{};

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

		return TotalValueCount;
	}


	// get the percentage of a kind of value in a field
	FloatType getPercentage(const FieldDataType value, std::size_t fieldidx = 0, bool useolap = true) const {
		std::size_t TotalCount{};
		std::size_t TotalValueCount{};

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

		return static_cast<FloatType>(TotalValueCount) / TotalCount;
	}


	// get number of a kind of flag in a flag field
	std::size_t getFlagCount(const FieldDataType flag, std::size_t fieldidx = 0, bool useolap = true) const {
		std::size_t TotalValueCount{};

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

		return TotalValueCount;
	}


	// get the percentage of a kind of flag value in a flag field
	FloatType getFlagPercentage(const FieldDataType flag, std::size_t fieldidx = 0, bool useolap = true) const {
		std::size_t TotalCount{};
		std::size_t TotalValueCount{};

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

		return static_cast<FloatType>(TotalValueCount) / TotalCount;
	}


	// print statistics of a flag field
	void printFlagStatistics(std::size_t fieldidx = 0, bool printpercentage = true, bool useolap = true) const {
		std::size_t TotalCount{};
		// get the number of different flags
		std::vector<std::pair<FieldDataType, std::size_t>> flagcount;

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
		// sort the vector
		std::sort(flagcount.begin(), flagcount.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
		// print info
		MPI_RANK(0);
		std::cout << "[Flag statistics]:" << std::endl;
		std::cout << "Flag | ";
		if (printpercentage) std::cout << "Percent% | ";
		std::cout << "Count" << std::endl;
		for (const auto& p : flagcount) {
			std::cout << std::setw(7) << std::left << static_cast<int>(p.first);
			if (printpercentage) std::cout << std::setw(11) << std::left << static_cast<FloatType>(p.second) / TotalCount * 100;
			std::cout << p.second << std::endl;
		}
		std::cout << std::endl;
	}


};