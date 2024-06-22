/* This file is part of FreeLB
 * 
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 * 
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

#pragma once

#include <stdexcept>

#include <omp.h>
#include <stdlib.h>
// #define Thread_Num 16
// // if define Thread_Nums in makefile
// // e.g., CXXFLAGS = -DThread_Nums=4 -fopenmp

// you can define Thread_Num in makefile:
// FLAGS += -DThread_Num=16 -fopenmp
#ifndef Thread_Num
// tells the compiler Thread_Num variable is defined somewhere else
// if used in code, need to define it
// you can define Thread_Num in .cpp file
// then use .ini file to store Thread_Num
extern int Thread_Num;
#endif


// float number type
// FLAGS        += -DFLOAT_TYPE=double
#ifdef FLOAT_TYPE
using FLOAT = FLOAT_TYPE;
#else
using FLOAT = double;
#endif


// users can decide whether to add the debug macros:
// FLAGS        += -D_FLB_DEBUG