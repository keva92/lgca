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

#ifndef LGCA_MODELS_H_
#define LGCA_MODELS_H_

#include "lgca_common.h"

namespace lgca {

// Struct for model-based values according to the number of lattice directions
template<Model model_>
struct ModelDescriptor;

// HPP model
template<>
struct ModelDescriptor<Model::HPP> {

    static constexpr unsigned int NUM_DIR = 4;

    // Inverse direction indices for each lattice direction
    static constexpr char INV_DIR       [NUM_DIR] = {   2,    3,    0,    1};

    // Mirrored direction indices for each lattice direction with respect to the x and y axis
    static constexpr char MIR_DIR_X     [NUM_DIR] = {   0,    3,    2,    1};
    static constexpr char MIR_DIR_Y     [NUM_DIR] = {   2,    1,    0,    3};

    // Lattice vector components in the different directions
    static constexpr Real LATTICE_VEC_X [NUM_DIR] = { 1.0,  0.0, -1.0,  0.0}; // = cos(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))
    static constexpr Real LATTICE_VEC_Y [NUM_DIR] = { 0.0,  1.0,  0.0, -1.0}; // = sin(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))

    // Memory offset to neighbor cells in the different directions for the propagation step
    int offset_to_neighbor_even         [NUM_DIR];
    int offset_to_neighbor_odd          [NUM_DIR];

    // Memory offset to related cells of the opposite boundary in the different directions in case
    // of periodic boundaries
    int offset_to_eastern_boundary_even [NUM_DIR];
    int offset_to_eastern_boundary_odd  [NUM_DIR];
    int offset_to_northern_boundary_even[NUM_DIR];
    int offset_to_northern_boundary_odd [NUM_DIR];
    int offset_to_western_boundary_even [NUM_DIR];
    int offset_to_western_boundary_odd  [NUM_DIR];
    int offset_to_southern_boundary_even[NUM_DIR];
    int offset_to_southern_boundary_odd [NUM_DIR];

    ModelDescriptor(const unsigned int dim_x, const unsigned int dim_y)
    {
        // Cell located in a row with even index value
        offset_to_neighbor_even         [0] = 1;
        offset_to_neighbor_even         [1] = dim_x;
        offset_to_neighbor_even         [2] = -1;
        offset_to_neighbor_even         [3] = -dim_x;

        offset_to_eastern_boundary_even [0] = 0;
        offset_to_eastern_boundary_even [1] = 0;
        offset_to_eastern_boundary_even [2] = dim_x;
        offset_to_eastern_boundary_even [3] = 0;

        offset_to_northern_boundary_even[0] = 0;
        offset_to_northern_boundary_even[1] = 0;
        offset_to_northern_boundary_even[2] = 0;
        offset_to_northern_boundary_even[3] = dim_x * dim_y;

        offset_to_western_boundary_even [0] = -dim_x;
        offset_to_western_boundary_even [1] = 0;
        offset_to_western_boundary_even [2] = 0;
        offset_to_western_boundary_even [3] = 0;

        offset_to_southern_boundary_even[0] = 0;
        offset_to_southern_boundary_even[1] = -dim_x * dim_y;
        offset_to_southern_boundary_even[2] = 0;
        offset_to_southern_boundary_even[3] = 0;

        // Cell located in a row with odd index value
        offset_to_neighbor_odd          [0] = offset_to_neighbor_even           [0];
        offset_to_neighbor_odd          [1] = offset_to_neighbor_even           [1];
        offset_to_neighbor_odd          [2] = offset_to_neighbor_even           [2];
        offset_to_neighbor_odd          [3] = offset_to_neighbor_even           [3];

        offset_to_eastern_boundary_odd  [0] = offset_to_eastern_boundary_even   [0];
        offset_to_eastern_boundary_odd  [1] = offset_to_eastern_boundary_even   [1];
        offset_to_eastern_boundary_odd  [2] = offset_to_eastern_boundary_even   [2];
        offset_to_eastern_boundary_odd  [3] = offset_to_eastern_boundary_even   [3];

        offset_to_northern_boundary_odd [0] = offset_to_northern_boundary_even  [0];
        offset_to_northern_boundary_odd [1] = offset_to_northern_boundary_even  [1];
        offset_to_northern_boundary_odd [2] = offset_to_northern_boundary_even  [2];
        offset_to_northern_boundary_odd [3] = offset_to_northern_boundary_even  [3];

        offset_to_western_boundary_odd  [0] = offset_to_western_boundary_even   [0];
        offset_to_western_boundary_odd  [1] = offset_to_western_boundary_even   [1];
        offset_to_western_boundary_odd  [2] = offset_to_western_boundary_even   [2];
        offset_to_western_boundary_odd  [3] = offset_to_western_boundary_even   [3];

        offset_to_southern_boundary_odd [0] = offset_to_southern_boundary_even  [0];
        offset_to_southern_boundary_odd [1] = offset_to_southern_boundary_even  [1];
        offset_to_southern_boundary_odd [2] = offset_to_southern_boundary_even  [2];
        offset_to_southern_boundary_odd [3] = offset_to_southern_boundary_even  [3];
    }
};

// FHP model
template<>
struct ModelDescriptor<Model::FHP> {

    static constexpr unsigned int NUM_DIR = 6;

    static constexpr Real SIN = sin(M_PI/3);

    // Inverse direction indices for each lattice direction
    static constexpr char INV_DIR       [NUM_DIR] = {   3,    4,    5,    0,    1,    2};

    // Mirrored direction indices for each lattice direction with respect to the x and y axis
    static constexpr char MIR_DIR_X     [NUM_DIR] = {   0,    5,    4,    3,    2,    1};
    static constexpr char MIR_DIR_Y     [NUM_DIR] = {   3,    2,    1,    0,    5,    4};

    // Lattice vector components in the different directions
    static constexpr Real LATTICE_VEC_X [NUM_DIR] = { 1.0,  0.5, -0.5, -1.0, -0.5,  0.5}; // = cos(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))
    static constexpr Real LATTICE_VEC_Y [NUM_DIR] = { 0.0,  SIN,  SIN,  0.0, -SIN, -SIN}; // = sin(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))

    // Memory offset to neighbor cells in the different directions for the propagation step
    // Note that for the FHP model there is a difference in the offsets depending on weather the
    // cell is located in a row with even or odd index.
    int offset_to_neighbor_even         [NUM_DIR];
    int offset_to_neighbor_odd          [NUM_DIR];

    // Memory offset to related cells of the opposite boundary in the different directions in case
    // of periodic boundaries
    int offset_to_eastern_boundary_even [NUM_DIR];
    int offset_to_eastern_boundary_odd  [NUM_DIR];
    int offset_to_northern_boundary_even[NUM_DIR];
    int offset_to_northern_boundary_odd [NUM_DIR];
    int offset_to_western_boundary_even [NUM_DIR];
    int offset_to_western_boundary_odd  [NUM_DIR];
    int offset_to_southern_boundary_even[NUM_DIR];
    int offset_to_southern_boundary_odd [NUM_DIR];

    ModelDescriptor(const unsigned int dim_x, const unsigned int dim_y)
    {
        // Cell located in a row with even index value
        offset_to_neighbor_even         [0] = 1;
        offset_to_neighbor_even         [1] = dim_x;
        offset_to_neighbor_even         [2] = dim_x - 1;
        offset_to_neighbor_even         [3] = -1;
        offset_to_neighbor_even         [4] = -dim_x - 1;
        offset_to_neighbor_even         [5] = -dim_x;

        offset_to_eastern_boundary_even [0] = 0;
        offset_to_eastern_boundary_even [1] = 0;
        offset_to_eastern_boundary_even [2] = dim_x;
        offset_to_eastern_boundary_even [3] = dim_x;
        offset_to_eastern_boundary_even [4] = dim_x;
        offset_to_eastern_boundary_even [5] = 0;

        offset_to_northern_boundary_even[0] = 0;
        offset_to_northern_boundary_even[1] = 0;
        offset_to_northern_boundary_even[2] = 0;
        offset_to_northern_boundary_even[3] = 0;
        offset_to_northern_boundary_even[4] = dim_x * dim_y;
        offset_to_northern_boundary_even[5] = dim_x * dim_y;

        offset_to_western_boundary_even [0] = -dim_x;
        offset_to_western_boundary_even [1] = 0;
        offset_to_western_boundary_even [2] = 0;
        offset_to_western_boundary_even [3] = 0;
        offset_to_western_boundary_even [4] = 0;
        offset_to_western_boundary_even [5] = 0;

        offset_to_southern_boundary_even[0] = 0;
        offset_to_southern_boundary_even[1] = -dim_x * dim_y;
        offset_to_southern_boundary_even[2] = -dim_x * dim_y + 1;
        offset_to_southern_boundary_even[3] = 0;
        offset_to_southern_boundary_even[4] = 0;
        offset_to_southern_boundary_even[5] = 0;

        // Cell located in a row with odd index value
        offset_to_neighbor_odd          [0] = 1;
        offset_to_neighbor_odd          [1] = dim_x + 1;
        offset_to_neighbor_odd          [2] = dim_x;
        offset_to_neighbor_odd          [3] = -1;
        offset_to_neighbor_odd          [4] = -dim_x;
        offset_to_neighbor_odd          [5] = -dim_x + 1;

        offset_to_eastern_boundary_odd  [0] = 0;
        offset_to_eastern_boundary_odd  [1] = 0;
        offset_to_eastern_boundary_odd  [2] = 0;
        offset_to_eastern_boundary_odd  [3] = dim_x;
        offset_to_eastern_boundary_odd  [4] = 0;
        offset_to_eastern_boundary_odd  [5] = 0;

        offset_to_northern_boundary_odd [0] = 0;
        offset_to_northern_boundary_odd [1] = 0;
        offset_to_northern_boundary_odd [2] = 0;
        offset_to_northern_boundary_odd [3] = 0;
        offset_to_northern_boundary_odd [4] = dim_x * dim_y;
        offset_to_northern_boundary_odd [5] = dim_x * dim_y;

        offset_to_western_boundary_odd  [0] = -dim_x;
        offset_to_western_boundary_odd  [1] = -dim_x;
        offset_to_western_boundary_odd  [2] = 0;
        offset_to_western_boundary_odd  [3] = 0;
        offset_to_western_boundary_odd  [4] = 0;
        offset_to_western_boundary_odd  [5] = -dim_x;

        offset_to_southern_boundary_odd [0] = 0;
        offset_to_southern_boundary_odd [1] = -dim_x * dim_y;
        offset_to_southern_boundary_odd [2] = -dim_x * dim_y;
        offset_to_southern_boundary_odd [3] = 0;
        offset_to_southern_boundary_odd [4] = 0;
        offset_to_southern_boundary_odd [5] = 0;
    }
};

} // namespace lgca

#endif /* LGCA_MODELS_H_ */
