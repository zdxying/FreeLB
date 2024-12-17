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

// block_reader.h
// read basicblock info from file and construct the objects
// parameters to construct BasicBlock

// T voxsize, const AABB<T, D>& aabb, const AABB<int, D>& idxblock, int blockid = -1
// or
// std::uint8_t level, T newvoxsize, int blockid, const AABB<T, D>& aabb, const AABB<int, D>& idxblock, const Vector<int, D>& mesh

// the file should be in the format:
/* # this is a comment
 * Refine: 0 / 1  # 0 for plain block, 1 for refined block
 * VoxelSize: 1.  # float or integer type
 * [Block]
 * Level: 0   # integer
 * AABB: {{minx, miny}, {maxx, maxy}}  # float or integer type, 2D or 3D
 * Mesh: {Nx, Ny}  # integer type, 2D or 3D
 * # empty line
 * [Block]
 * # next block
*/

#pragma once

#include <fstream>
#include <stdio.h>
#include <string>
#include <cstdint>

#include "geometry/basic_geometry.h"
#include "utils/util.h"


template <typename T, unsigned int D>
class Vector;


template <typename T>
T StrToType(const std::string& str) {
	// trim spaces
	std::string substr = util::trim(str);
	// convert to T
	return T(std::stod(substr));
}

template <typename T>
T StrToTypeInt(const std::string& str) {
	// trim spaces
	std::string substr = util::trim(str);
	// convert to int
	return std::stoi(substr);
}

template <typename T, unsigned int D>
Vector<T,D> StrToVector(const std::string& str) {
	std::string sub_str = str.substr(str.find("{") + 1, str.find("}") - str.find("{") - 1);
	T x = StrToType<T>(sub_str.substr(0, sub_str.find(",")));
	sub_str = sub_str.substr(sub_str.find(",") + 1);
	if constexpr (D == 2) {
		T y = StrToType<T>(sub_str.substr(0, sub_str.find("}")));	
		return Vector<T, 2>(x, y);
	} else if constexpr (D == 3) {
		T y = StrToType<T>(sub_str.substr(0, sub_str.find(",")));	
		sub_str = sub_str.substr(sub_str.find(",") + 1);
		T z = StrToType<T>(sub_str.substr(0, sub_str.find("}")));
		return Vector<T, 3>(x, y, z);
	}
}

template <typename T, unsigned int D>
Vector<T,D> StrToVectorInt(const std::string& str) {
	std::string sub_str = str.substr(str.find("{") + 1, str.find("}") - str.find("{") - 1);
	T x = StrToTypeInt<T>(sub_str.substr(0, sub_str.find(",")));
	sub_str = sub_str.substr(sub_str.find(",") + 1);
	if constexpr (D == 2) {
		T y = StrToTypeInt<T>(sub_str.substr(0, sub_str.find("}")));	
		return Vector<T, 2>(x, y);
	} else if constexpr (D == 3) {
		T y = StrToTypeInt<T>(sub_str.substr(0, sub_str.find(",")));	
		sub_str = sub_str.substr(sub_str.find(",") + 1);
		T z = StrToTypeInt<T>(sub_str.substr(0, sub_str.find("}")));
		return Vector<T, 3>(x, y, z);
	}
}


// read multiple 2D basicblocks' info from file
template <typename T>
class BlockReader2D {
	private:
	std::vector<BasicBlock<T, 2>> blocks;
	T VoxelSize;
	// an aabb including all blocks
	BasicBlock<T, 2> BaseBlock;

	std::uint8_t MaxLevel;

	void getBaseBlock(){
		T minx = std::numeric_limits<T>::max();
		T miny = std::numeric_limits<T>::max();
		T maxx = std::numeric_limits<T>::min();
		T maxy = std::numeric_limits<T>::min();
		int idxminx{};
		int idxminy{};
		int idxmaxx{};
		int idxmaxy{};
		for (const BasicBlock<T, 2>& block : blocks) {
			if (block.getMin()[0] < minx) minx = block.getMin()[0];
			if (block.getMin()[1] < miny) miny = block.getMin()[1];
			if (block.getMax()[0] > maxx) maxx = block.getMax()[0];
			if (block.getMax()[1] > maxy) maxy = block.getMax()[1];

			if (block.getIdxBlock().getMin()[0] < idxminx) idxminx = block.getIdxBlock().getMin()[0];
			if (block.getIdxBlock().getMin()[1] < idxminy) idxminy = block.getIdxBlock().getMin()[1];
			if (block.getIdxBlock().getMax()[0] > idxmaxx) idxmaxx = block.getIdxBlock().getMax()[0];
			if (block.getIdxBlock().getMax()[1] > idxmaxy) idxmaxy = block.getIdxBlock().getMax()[1];
		}
		AABB<T, 2> aabb{Vector<T, 2>{minx, miny}, Vector<T, 2>{maxx, maxy}};
		AABB<int, 2> idxaabb{Vector<int, 2>{idxminx, idxminy}, Vector<int, 2>{idxmaxx, idxmaxy}};
		BaseBlock = BasicBlock<T, 2>(VoxelSize, aabb, idxaabb, -1);
	}

	void getMaxLevel(){
		for (const BasicBlock<T, 2>& block : blocks) {
			if (block.getLevel() > MaxLevel) MaxLevel = block.getLevel();
		}
	}

	public:
	BlockReader2D(const std::string& filename) : VoxelSize(), MaxLevel(std::uint8_t{}), BaseBlock() {
		std::ifstream file(filename);
		if (!file) {
			std::cerr << "Error: cannot open file " << filename << std::endl;
			exit(1);
		}
		std::string line, section;

		int blockid = -1;
		bool refine = false;
		std::uint8_t level{};
		{};
		AABB<T, 2> aabb{};
		Vector<int, 2> mesh{};

		// read Refine
		while (std::getline(file, line)) {
			// ignore empty lines
			if (line.empty()) continue;
			// ignore comments
			if (line[0] == '#') continue;

			if (line.find("Refine") != std::string::npos) {
				std::string refine_str = line.substr(line.find(":") + 1);
				// trim spaces
				refine_str = util::trim(refine_str);
				// convert to bool
				refine = (refine_str == "1");
				break;
			}
		}
		// read VoxelSize
		while (std::getline(file, line)) {
			// ignore empty lines
			if (line.empty()) continue;
			// ignore comments
			if (line[0] == '#') continue;

			if (line.find("VoxelSize") != std::string::npos) {
				VoxelSize = StrToType<T>(line.substr(line.find(":") + 1));
				break;
			}
		}

		while (std::getline(file, line)) {
			if (line.empty()) continue;
			if (line[0] == '#') continue;
			if (line[0] == '[') {
				// read a new block
				++blockid;
				level = std::uint8_t{};
				aabb = AABB<T, 2>{};
				mesh = Vector<int, 2>{};
				continue;
			} else {
				// read block info
				// read Level
				if (refine){
					if (line.find("Level") != std::string::npos) {
						level = StrToTypeInt<std::uint8_t>(line.substr(line.find(":") + 1));
					}
				}
				// read AABB: {{minx, miny}, {maxx, maxy}}
				if (line.find("AABB") != std::string::npos) {
					std::string sub_str = line.substr(line.find("}") + 1);

					Vector<T, 2> min = StrToVector<T, 2>(line.substr(line.find("{") + 1));
					Vector<T, 2> max = StrToVector<T, 2>(sub_str.substr(sub_str.find("{")));
					// convert to T
					aabb = AABB<T, 2>(min, max);
				}
				// read Mesh: {Nx, Ny}
				if (line.find("Mesh") != std::string::npos) {
					mesh = StrToVectorInt<int, 2>(line.substr(line.find(":") + 1));
					// index block
					AABB<int, 2> idxaabb{Vector<int, 2>{}, mesh - Vector<int, 2>{1}};
					// construct BasicBlock
					if (refine) {
						int coeff = int(std::pow(2, int(level)));
						Vector<int, 2> refmesh = mesh * coeff;				
						blocks.emplace_back(level, VoxelSize, blockid, aabb, idxaabb, refmesh);
					} else {
						blocks.emplace_back(VoxelSize, aabb, idxaabb, blockid);
					}
				}
			}
		}
		
		getMaxLevel();
		getBaseBlock();
	}

	const std::vector<BasicBlock<T, 2>>& getBlocks() const {
		return blocks;
	}
	const BasicBlock<T, 2>& getBaseBlock() const {
		return BaseBlock;
	}
	BasicBlock<T, 2> getBasicBlock() const {
		return BaseBlock.getExtBlock(1);
	}
	std::uint8_t getMaxLevel() const {
		return MaxLevel;
	}

};

// read multiple 3D basicblocks' info from file
template <typename T>
class BlockReader3D {
	private:
	T VoxelSize;
	std::uint8_t MaxLevel;
	// an aabb including all blocks
	BasicBlock<T, 3> BaseBlock;

	std::vector<BasicBlock<T, 3>> blocks;
	std::vector<int> overlaps;

	void getBaseBlock(){
		T minx = std::numeric_limits<T>::max();
		T miny = std::numeric_limits<T>::max();
		T minz = std::numeric_limits<T>::max();
		T maxx = std::numeric_limits<T>::min();
		T maxy = std::numeric_limits<T>::min();
		T maxz = std::numeric_limits<T>::min();
		int idxminx{};
		int idxminy{};
		int idxminz{};
		int idxmaxx{};
		int idxmaxy{};
		int idxmaxz{};
		for (const BasicBlock<T, 3>& block : blocks) {
			if (block.getMin()[0] < minx) minx = block.getMin()[0];
			if (block.getMin()[1] < miny) miny = block.getMin()[1];
			if (block.getMin()[2] < minz) minz = block.getMin()[2];
			if (block.getMax()[0] > maxx) maxx = block.getMax()[0];
			if (block.getMax()[1] > maxy) maxy = block.getMax()[1];
			if (block.getMax()[2] > maxz) maxz = block.getMax()[2];

			if (block.getIdxBlock().getMin()[0] < idxminx) idxminx = block.getIdxBlock().getMin()[0];
			if (block.getIdxBlock().getMin()[1] < idxminy) idxminy = block.getIdxBlock().getMin()[1];
			if (block.getIdxBlock().getMin()[2] < idxminz) idxminz = block.getIdxBlock().getMin()[2];
			if (block.getIdxBlock().getMax()[0] > idxmaxx) idxmaxx = block.getIdxBlock().getMax()[0];
			if (block.getIdxBlock().getMax()[1] > idxmaxy) idxmaxy = block.getIdxBlock().getMax()[1];
			if (block.getIdxBlock().getMax()[2] > idxmaxz) idxmaxz = block.getIdxBlock().getMax()[2];
		}
		AABB<T, 3> aabb{Vector<T, 3>{minx, miny, minz}, Vector<T, 3>{maxx, maxy, maxz}};
		AABB<int, 3> idxaabb{Vector<int, 3>{idxminx, idxminy, idxminz}, Vector<int, 3>{idxmaxx, idxmaxy, idxmaxz}};
		BaseBlock = BasicBlock<T, 3>(VoxelSize, aabb, idxaabb, -1);
	}

	void getMaxLevel(){
		for (const BasicBlock<T, 3>& block : blocks) {
			if (block.getLevel() > MaxLevel) MaxLevel = block.getLevel();
		}
	}

	public:
	BlockReader3D(const std::string& filename) : VoxelSize(), MaxLevel(std::uint8_t{}), BaseBlock() {
		std::ifstream file(filename);
		if (!file) {
			std::cerr << "Error: cannot open file " << filename << std::endl;
			exit(1);
		}
		std::string line, section;

		int blockid = -1;
		bool refine = false;
		std::uint8_t level{};
		{};
		AABB<T, 3> aabb{};
		Vector<int, 3> mesh{};

		// read Refine
		while (std::getline(file, line)) {
			// ignore empty lines
			if (line.empty()) continue;
			// ignore comments
			if (line[0] == '#') continue;

			if (line.find("Refine") != std::string::npos) {
				std::string refine_str = line.substr(line.find(":") + 1);
				// trim spaces
				refine_str = util::trim(refine_str);
				// convert to bool
				refine = (refine_str == "1");
				break;
			}
		}
		// read VoxelSize
		while (std::getline(file, line)) {
			// ignore empty lines
			if (line.empty()) continue;
			// ignore comments
			if (line[0] == '#') continue;

			if (line.find("VoxelSize") != std::string::npos) {
				VoxelSize = StrToType<T>(line.substr(line.find(":") + 1));
				break;
			}
		}

		while (std::getline(file, line)) {
			if (line.empty()) continue;
			if (line[0] == '#') continue;
			if (line[0] == '[') {
				// read a new block
				++blockid;
				level = std::uint8_t{};
				aabb = AABB<T, 3>{};
				mesh = Vector<int, 3>{};
				continue;
			} else {
				// read block info
				// read overlaps
				if (line.find("Overlap") != std::string::npos) {
					overlaps.push_back(StrToTypeInt<int>(line.substr(line.find(":") + 1)));
				}
				// read Level
				if (refine){
					if (line.find("Level") != std::string::npos) {
						level = StrToTypeInt<std::uint8_t>(line.substr(line.find(":") + 1));
					}
				}
				// read AABB: {{minx, miny}, {maxx, maxy}}
				if (line.find("AABB") != std::string::npos) {
					std::string sub_str = line.substr(line.find("}") + 1);

					Vector<T, 3> min = StrToVector<T, 3>(line.substr(line.find("{") + 1));
					Vector<T, 3> max = StrToVector<T, 3>(sub_str.substr(sub_str.find("{")));
					// convert to T
					aabb = AABB<T, 3>(min, max);
				}
				// read Mesh: {Nx, Ny}
				if (line.find("Mesh") != std::string::npos) {
					mesh = StrToVectorInt<int, 3>(line.substr(line.find(":") + 1));
					// index block
					AABB<int, 3> idxaabb{Vector<int, 3>{}, mesh - Vector<int, 3>{1}};
					// construct BasicBlock
					if (refine) {
						int coeff = int(std::pow(2, int(level)));
						Vector<int, 3> refmesh = mesh * coeff;				
						blocks.emplace_back(level, VoxelSize, blockid, aabb, idxaabb, refmesh);
					} else {
						blocks.emplace_back(VoxelSize, aabb, idxaabb, blockid);
					}
				}
			}
		}
		
		getMaxLevel();
		getBaseBlock();
	}

	const std::vector<BasicBlock<T, 3>>& getBlocks() const {
		return blocks;
	}
	const std::vector<int>& getOverlaps() const {
		return overlaps;
	}
	const BasicBlock<T, 3>& getBaseBlock() const {
		return BaseBlock;
	}
	BasicBlock<T, 3> getBasicBlock() const {
		return BaseBlock.getExtBlock(1);
	}
	std::uint8_t getMaxLevel() const {
		return MaxLevel;
	}

};


