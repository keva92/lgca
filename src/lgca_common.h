/*
 * This file is part of LGCA, an implementation of a Lattice Gas Cellular Automaton
 * (https://github.com/keva92/lgca).
 *
 * Copyright (c) 2015-2017 Kerstin Vater, Niklas Kühl, Christian F. Janßen.
 *
 * LGCA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * LGCA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with lgca. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LGCA_COMMON_H_
#define LGCA_COMMON_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cassert>
#include <chrono>

namespace lgca {

#define WORD_SIZE 8;

using std::cout;
using std::endl;
using std::flush;
using std::string;

using std::chrono::steady_clock;
using std::chrono::duration;

// Define floating-point precision
typedef float Real;

// Number of lattice directions defines the different lattice gas models
enum class Model {
    HPP = 4,
    FHP = 6
};

enum class CellType {
    FLUID         = 0,
    SOLID_NO_SLIP = 1,
    SOLID_SLIP    = 2
};

} // namespace lgca

#endif /* LGCA_COMMON_H_ */
