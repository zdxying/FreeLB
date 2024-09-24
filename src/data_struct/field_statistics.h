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
	FieldDataType getAverage(std::size_t fieldidx = 0) const {
		FieldDataType TotalSum{};
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(+ : TotalSum)
#endif
		for (const auto& blockfield : BFM.getBlockFields()) {
			const auto& field = blockfield.getFieldType().getField(fieldidx);
			FieldDataType sum{};
			const std::size_t size = field.size();
			for (std::size_t i = 0; i < size; ++i) {
				sum += field[i];
			}
			sum /= size;
			TotalSum += sum;
		}
		return TotalSum / BFM.size();
	}

	// get the maximum value of a field
	FieldDataType getMax(std::size_t fieldidx = 0) const {
		FieldDataType TotalMax{};
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : TotalMax)
#endif
		for (const auto& blockfield : BFM.getBlockFields()) {
			const auto& field = blockfield.getFieldType().getField(fieldidx);
			FieldDataType max{};
			const std::size_t size = field.size();
			for (std::size_t i = 0; i < size; ++i) {
				max = max > field[i] ? max : field[i];
			}
			TotalMax = TotalMax > max ? TotalMax : max;
		}
		return TotalMax;
	}

	// get the minimum value of a field
	FieldDataType getMin(std::size_t fieldidx = 0) const {
		FieldDataType TotalMin{};
#ifdef _OPENMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(min : TotalMin)
#endif
		for (const auto& blockfield : BFM.getBlockFields()) {
			const auto& field = blockfield.getFieldType().getField(fieldidx);
			FieldDataType min{};
			const std::size_t size = field.size();
			for (std::size_t i = 0; i < size; ++i) {
				min = min < field[i] ? min : field[i];
			}
			TotalMin = TotalMin < min ? TotalMin : min;
		}
		return TotalMin;
	}

	// get the percentage of a kind of value in a field
	FloatType getPercentage(FieldDataType value, std::size_t fieldidx = 0) const {
		std::size_t TotalCount{};
		std::size_t TotalValueCount{};
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
		return static_cast<FloatType>(TotalValueCount) / TotalCount;
	}

};