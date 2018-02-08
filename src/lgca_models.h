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

    // Four directions
    //     1
    //     |
    // 2 -   - 0
    //     |
    //     3
    static constexpr unsigned int NUM_DIR = 4;

    // Inverse direction indices for each lattice direction
    static constexpr char INV_DIR       [NUM_DIR] = {   2,    3,    0,    1};

    // Mirrored direction indices for each lattice direction with respect to the x and y axis
    static constexpr char MIR_DIR_X     [NUM_DIR] = {   0,    3,    2,    1};
    static constexpr char MIR_DIR_Y     [NUM_DIR] = {   2,    1,    0,    3};

    // Lattice vector components in the different directions
    static constexpr Real LATTICE_VEC_X [NUM_DIR] = { 1.0,  0.0, -1.0,  0.0}; // = cos(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))
    static constexpr Real LATTICE_VEC_Y [NUM_DIR] = { 0.0,  1.0,  0.0, -1.0}; // = sin(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))

    // TODO Collision table
    static constexpr unsigned char COLLISION_LUT[1 << NUM_DIR] = { };

    // TODO Bounce back table
    static constexpr unsigned char BB_LUT[1 << NUM_DIR] = { };

    // TODO Bounce forward tables
    static constexpr unsigned char BF_X_LUT[1 << NUM_DIR] = { };
    static constexpr unsigned char BF_Y_LUT[1 << NUM_DIR] = { };

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

    static inline void collide(unsigned char* node_state_in, unsigned char* node_state_out, const bool p)
    {
//        node_state_out[0] = COLLISION_LUT[node_state_in[0]];

        node_state_out[0] = node_state_in[0]
                - (node_state_in[0] * node_state_in[2] * (1 - node_state_in[1]) * (1 - node_state_in[3]))
                + (node_state_in[1] * node_state_in[3] * (1 - node_state_in[0]) * (1 - node_state_in[2]));

        node_state_out[1] = node_state_in[1]
                - (node_state_in[1] * node_state_in[3] * (1 - node_state_in[0]) * (1 - node_state_in[2]))
                + (node_state_in[0] * node_state_in[2] * (1 - node_state_in[1]) * (1 - node_state_in[3]));

        node_state_out[2] = node_state_in[2]
                - (node_state_in[0] * node_state_in[2] * (1 - node_state_in[1]) * (1 - node_state_in[3]))
                + (node_state_in[1] * node_state_in[3] * (1 - node_state_in[0]) * (1 - node_state_in[2]));

        node_state_out[3] = node_state_in[3]
                - (node_state_in[1] * node_state_in[3] * (1 - node_state_in[0]) * (1 - node_state_in[2]))
                + (node_state_in[0] * node_state_in[2] * (1 - node_state_in[1]) * (1 - node_state_in[3]));

//        // Collision case 1.
//        if ((node_state_in[0] == 0) &&
//            (node_state_in[1] == 1) &&
//            (node_state_in[2] == 0) &&
//            (node_state_in[3] == 1)) {

//            node_state_out[0] = 1;
//            node_state_out[1] = 0;
//            node_state_out[2] = 1;
//            node_state_out[3] = 0;

//            return;
//        }

//        // Collision case 2.
//        if ((node_state_in[0] == 1) &&
//            (node_state_in[1] == 0) &&
//            (node_state_in[2] == 1) &&

//            (node_state_in[3] == 0)) {

//            node_state_out[0] = 0;
//            node_state_out[1] = 1;
//            node_state_out[2] = 0;
//            node_state_out[3] = 1;

//            return;
//        }
    }

    static inline void bounce_back(unsigned char* node_state_in, unsigned char* node_state_out)
    {
//        node_state_out[0] = BB_LUT[node_state_in[0]];

#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

            // Exchange the states of the nodes with the the states of the inverse directions
            node_state_out[dir] = node_state_in[INV_DIR[dir]];
        }
    }

    static inline void bounce_forward_x(unsigned char* node_state_in, unsigned char* node_state_out)
    {
//        node_state_out[0] = BF_X_LUT[node_state_in[0]];

#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

                // Exchange the states of the nodes with the the states of the mirrored
                // directions along the x axis
                node_state_out[dir] = node_state_in[MIR_DIR_X[dir]];
        }
    }

    static inline void bounce_forward_y(unsigned char* node_state_in, unsigned char* node_state_out)
    {
//        node_state_out[0] = BF_Y_LUT[node_state_in[0]];

#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

                // Exchange the states of the nodes with the the states of
                // the mirrored directions along the y axis
                node_state_out[dir] = node_state_in[MIR_DIR_Y[dir]];
        }
    }

}; // struct ModelDescriptor<Model::HPP>

// FHP-I model
template<>
struct ModelDescriptor<Model::FHP_I> {

    // Six directions
    //   2    1
    //    \  /
    // 3 -   - 0
    //    / \
    //   4   5
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

    // Collision table
    static constexpr unsigned char COLLISION_LUT[1 << NUM_DIR] = {
            0,  1,  2,  3,  4,  5,  6,  7,
            8, 18, 10, 11, 12, 13, 14, 15,
           16, 17, 36, 19, 20, 42, 22, 23,
           24, 25, 26, 27, 28, 29, 30, 31,
           32, 33, 34, 35,  9, 37, 38, 39,
           40, 41, 21, 43, 44, 45, 46, 47,
           48, 49, 50, 51, 52, 53, 54, 55,
           56, 57, 58, 59, 60, 61, 61, 63};

    // Bounce back table
    static constexpr unsigned char BB_LUT[1 << NUM_DIR] = {
            0,  8, 16, 24, 32, 40, 48, 56,
            1,  9, 17, 25, 33, 41, 49, 57,
            2, 10, 18, 26, 34, 42, 50, 58,
            3, 11, 19, 27, 35, 43, 51, 59,
            4, 12, 20, 28, 36, 44, 52, 60,
            5, 13, 21, 29, 37, 45, 53, 61,
            6, 14, 22, 30, 38, 46, 54, 62,
            7, 15, 23, 31, 39, 47, 55, 63};

    // TODO Bounce forward tables
    static constexpr unsigned char BF_X_LUT[1 << NUM_DIR] = { };
    static constexpr unsigned char BF_Y_LUT[1 << NUM_DIR] = { };

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

    static inline void collide(unsigned char* node_state_in, unsigned char* node_state_out, const bool p)
    {
        unsigned char a = node_state_in[1];
        unsigned char b = node_state_in[2];
        unsigned char c = node_state_in[3];
        unsigned char d = node_state_in[4];
        unsigned char e = node_state_in[5];
        unsigned char f = node_state_in[0];

        unsigned char db1 = (a&d&~(b|c|e|f));
        unsigned char db2 = (b&e&~(a|c|d|f));
        unsigned char db3 = (c&f&~(a|b|d|e));

        unsigned char triple = (a^b)&(b^c)&(c^d)&(d^e)&(e^f);

        unsigned char cha = ((triple|db1|(p&db2)|((~p)&db3))&(~0));
        unsigned char chd = ((triple|db1|(p&db2)|((~p)&db3))&(~0));
        unsigned char chb = ((triple|db2|(p&db3)|((~p)&db1))&(~0));
        unsigned char che = ((triple|db2|(p&db3)|((~p)&db1))&(~0));
        unsigned char chc = ((triple|db3|(p&db1)|((~p)&db2))&(~0));
        unsigned char chf = ((triple|db3|(p&db1)|((~p)&db2))&(~0));

        node_state_out[1] = node_state_in[1]^cha;
        node_state_out[2] = node_state_in[2]^chb;
        node_state_out[3] = node_state_in[3]^chc;
        node_state_out[4] = node_state_in[4]^chd;
        node_state_out[5] = node_state_in[5]^che;
        node_state_out[0] = node_state_in[0]^chf;
    }

    static inline void bounce_back(unsigned char* node_state_in, unsigned char* node_state_out)
    {
#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

            // Exchange the states of the nodes with the the states of the inverse directions
            node_state_out[dir] = node_state_in[INV_DIR[dir]];
        }
    }

    static inline void bounce_forward_x(unsigned char* node_state_in, unsigned char* node_state_out)
    {
#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

                // Exchange the states of the nodes with the the states of the mirrored
                // directions along the x axis
                node_state_out[dir] = node_state_in[MIR_DIR_X[dir]];
        }
    }

    static inline void bounce_forward_y(unsigned char* node_state_in, unsigned char* node_state_out)
    {
#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

                // Exchange the states of the nodes with the the states of
                // the mirrored directions along the y axis
                node_state_out[dir] = node_state_in[MIR_DIR_Y[dir]];
        }
    }

}; // struct ModelDescriptor<Model::FHP_I>

// FHP-II model
template<>
struct ModelDescriptor<Model::FHP_II> {

    // Six direction + rest particle
    //   2    1
    //    \  /
    // 3 - 6 - 0
    //    / \
    //   4   5
    static constexpr unsigned int NUM_DIR = 7;

    static constexpr Real SIN = sin(M_PI/3);

    // Inverse direction indices for each lattice direction
    static constexpr char INV_DIR       [NUM_DIR] = {   3,    4,    5,    0,    1,    2,    6};

    // Mirrored direction indices for each lattice direction with respect to the x and y axis
    static constexpr char MIR_DIR_X     [NUM_DIR] = {   0,    5,    4,    3,    2,    1,    6};
    static constexpr char MIR_DIR_Y     [NUM_DIR] = {   3,    2,    1,    0,    5,    4,    6};

    // Lattice vector components in the different directions
    static constexpr Real LATTICE_VEC_X [NUM_DIR] = { 1.0,  0.5, -0.5, -1.0, -0.5,  0.5,  0.0}; // = cos(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))
    static constexpr Real LATTICE_VEC_Y [NUM_DIR] = { 0.0,  SIN,  SIN,  0.0, -SIN, -SIN,  0.0}; // = sin(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))

    // TODO Collision table
    static constexpr unsigned char COLLISION_LUT[1 << NUM_DIR] = { };

    // TODO Bounce back table
    static constexpr unsigned char BB_LUT[1 << NUM_DIR] = { };

    // TODO Bounce forward tables
    static constexpr unsigned char BF_X_LUT[1 << NUM_DIR] = { };
    static constexpr unsigned char BF_Y_LUT[1 << NUM_DIR] = { };

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
        offset_to_neighbor_even         [6] = 0;

        offset_to_eastern_boundary_even [0] = 0;
        offset_to_eastern_boundary_even [1] = 0;
        offset_to_eastern_boundary_even [2] = dim_x;
        offset_to_eastern_boundary_even [3] = dim_x;
        offset_to_eastern_boundary_even [4] = dim_x;
        offset_to_eastern_boundary_even [5] = 0;
        offset_to_eastern_boundary_even [6] = 0;

        offset_to_northern_boundary_even[0] = 0;
        offset_to_northern_boundary_even[1] = 0;
        offset_to_northern_boundary_even[2] = 0;
        offset_to_northern_boundary_even[3] = 0;
        offset_to_northern_boundary_even[4] = dim_x * dim_y;
        offset_to_northern_boundary_even[5] = dim_x * dim_y;
        offset_to_northern_boundary_even[6] = 0;

        offset_to_western_boundary_even [0] = -dim_x;
        offset_to_western_boundary_even [1] = 0;
        offset_to_western_boundary_even [2] = 0;
        offset_to_western_boundary_even [3] = 0;
        offset_to_western_boundary_even [4] = 0;
        offset_to_western_boundary_even [5] = 0;
        offset_to_western_boundary_even [6] = 0;

        offset_to_southern_boundary_even[0] = 0;
        offset_to_southern_boundary_even[1] = -dim_x * dim_y;
        offset_to_southern_boundary_even[2] = -dim_x * dim_y + 1;
        offset_to_southern_boundary_even[3] = 0;
        offset_to_southern_boundary_even[4] = 0;
        offset_to_southern_boundary_even[5] = 0;
        offset_to_southern_boundary_even[6] = 0;

        // Cell located in a row with odd index value
        offset_to_neighbor_odd          [0] = 1;
        offset_to_neighbor_odd          [1] = dim_x + 1;
        offset_to_neighbor_odd          [2] = dim_x;
        offset_to_neighbor_odd          [3] = -1;
        offset_to_neighbor_odd          [4] = -dim_x;
        offset_to_neighbor_odd          [5] = -dim_x + 1;
        offset_to_neighbor_odd          [6] = 0;

        offset_to_eastern_boundary_odd  [0] = 0;
        offset_to_eastern_boundary_odd  [1] = 0;
        offset_to_eastern_boundary_odd  [2] = 0;
        offset_to_eastern_boundary_odd  [3] = dim_x;
        offset_to_eastern_boundary_odd  [4] = 0;
        offset_to_eastern_boundary_odd  [5] = 0;
        offset_to_eastern_boundary_odd  [6] = 0;

        offset_to_northern_boundary_odd [0] = 0;
        offset_to_northern_boundary_odd [1] = 0;
        offset_to_northern_boundary_odd [2] = 0;
        offset_to_northern_boundary_odd [3] = 0;
        offset_to_northern_boundary_odd [4] = dim_x * dim_y;
        offset_to_northern_boundary_odd [5] = dim_x * dim_y;
        offset_to_northern_boundary_odd [6] = 0;

        offset_to_western_boundary_odd  [0] = -dim_x;
        offset_to_western_boundary_odd  [1] = -dim_x;
        offset_to_western_boundary_odd  [2] = 0;
        offset_to_western_boundary_odd  [3] = 0;
        offset_to_western_boundary_odd  [4] = 0;
        offset_to_western_boundary_odd  [5] = -dim_x;
        offset_to_western_boundary_odd  [6] = 0;

        offset_to_southern_boundary_odd [0] = 0;
        offset_to_southern_boundary_odd [1] = -dim_x * dim_y;
        offset_to_southern_boundary_odd [2] = -dim_x * dim_y;
        offset_to_southern_boundary_odd [3] = 0;
        offset_to_southern_boundary_odd [4] = 0;
        offset_to_southern_boundary_odd [5] = 0;
        offset_to_southern_boundary_odd [6] = 0;
    }

    static inline void collide(unsigned char* node_state_in, unsigned char* node_state_out, const bool p)
    {
        unsigned char a = node_state_in[1];
        unsigned char b = node_state_in[2];
        unsigned char c = node_state_in[3];
        unsigned char d = node_state_in[4];
        unsigned char e = node_state_in[5];
        unsigned char f = node_state_in[0];
        unsigned char r = node_state_in[6];

        unsigned char db1 = (a&d&~(b|c|e|f));
        unsigned char db2 = (b&e&~(a|c|d|f));
        unsigned char db3 = (c&f&~(a|b|d|e));

        unsigned char triple = (a^b)&(b^c)&(c^d)&(d^e)&(e^f);

        unsigned char ra = (r&a&~(b|c|d|e|f));
        unsigned char rb = (r&b&~(a|c|d|e|f));
        unsigned char rc = (r&c&~(a|b|d|e|f));
        unsigned char rd = (r&d&~(a|b|c|e|f));
        unsigned char re = (r&e&~(a|b|c|d|f));
        unsigned char rf = (r&f&~(a|b|c|d|e));

        unsigned char ra2 = (f&b&~(r|a|c|d|e));
        unsigned char rb2 = (a&c&~(r|b|d|e|f));
        unsigned char rc2 = (b&d&~(r|a|c|e|f));
        unsigned char rd2 = (c&e&~(r|a|b|d|f));
        unsigned char re2 = (d&f&~(r|a|b|c|e));
        unsigned char rf2 = (e&a&~(r|b|c|d|f));

        unsigned char cha = ((triple|db1|(p&db2)|((~p)&db3)|ra|rb|rf|ra2|rb2|rf2)&(~0));
        unsigned char chd = ((triple|db1|(p&db2)|((~p)&db3)|rd|rc|re|rd2|rc2|re2)&(~0));
        unsigned char chb = ((triple|db2|(p&db3)|((~p)&db1)|rb|ra|rc|rb2|ra2|rc2)&(~0));
        unsigned char che = ((triple|db2|(p&db3)|((~p)&db1)|re|rd|rf|re2|rd2|rf2)&(~0));
        unsigned char chc = ((triple|db3|(p&db1)|((~p)&db2)|rc|rb|rd|rc2|rb2|rd2)&(~0));
        unsigned char chf = ((triple|db3|(p&db1)|((~p)&db2)|rf|ra|re|rf2|ra2|re2)&(~0));
        unsigned char chr = ((ra|rb|rc|rd|re|rf|ra2|rb2|rc2|rd2|re2|rf2)&(~0));

        node_state_out[1] = node_state_in[1]^cha;
        node_state_out[2] = node_state_in[2]^chb;
        node_state_out[3] = node_state_in[3]^chc;
        node_state_out[4] = node_state_in[4]^chd;
        node_state_out[5] = node_state_in[5]^che;
        node_state_out[0] = node_state_in[0]^chf;
        node_state_out[6] = node_state_in[6]^chr;
    }

    static inline void bounce_back(unsigned char* node_state_in, unsigned char* node_state_out)
    {
#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

            // Exchange the states of the nodes with the the states of the inverse directions
            node_state_out[dir] = node_state_in[INV_DIR[dir]];
        }
    }

    static inline void bounce_forward_x(unsigned char* node_state_in, unsigned char* node_state_out)
    {
#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

                // Exchange the states of the nodes with the the states of the mirrored
                // directions along the x axis
                node_state_out[dir] = node_state_in[MIR_DIR_X[dir]];
        }
    }

    static inline void bounce_forward_y(unsigned char* node_state_in, unsigned char* node_state_out)
    {
#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

                // Exchange the states of the nodes with the the states of
                // the mirrored directions along the y axis
                node_state_out[dir] = node_state_in[MIR_DIR_Y[dir]];
        }
    }

}; // struct ModelDescriptor<Model::FHP_II>

// FHP-III model
template<>
struct ModelDescriptor<Model::FHP_III> {

    // Six direction + rest particle
    //   2    1
    //    \  /
    // 3 - 6 - 0
    //    / \
    //   4   5
    static constexpr unsigned int NUM_DIR = 7;

    static constexpr Real SIN = sin(M_PI/3);

    // Inverse direction indices for each lattice direction
    static constexpr char INV_DIR       [NUM_DIR] = {   3,    4,    5,    0,    1,    2,    6};

    // Mirrored direction indices for each lattice direction with respect to the x and y axis
    static constexpr char MIR_DIR_X     [NUM_DIR] = {   0,    5,    4,    3,    2,    1,    6};
    static constexpr char MIR_DIR_Y     [NUM_DIR] = {   3,    2,    1,    0,    5,    4,    6};

    // Lattice vector components in the different directions
    static constexpr Real LATTICE_VEC_X [NUM_DIR] = { 1.0,  0.5, -0.5, -1.0, -0.5,  0.5,  0.0}; // = cos(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))
    static constexpr Real LATTICE_VEC_Y [NUM_DIR] = { 0.0,  SIN,  SIN,  0.0, -SIN, -SIN,  0.0}; // = sin(2.0 * M_PI / ((Real) num_dir_) * ((Real) dir))

    // TODO Collision table
    static constexpr unsigned char COLLISION_LUT[1 << NUM_DIR] = { };

    // TODO Bounce back table
    static constexpr unsigned char BB_LUT[1 << NUM_DIR] = { };

    // TODO Bounce forward tables
    static constexpr unsigned char BF_X_LUT[1 << NUM_DIR] = { };
    static constexpr unsigned char BF_Y_LUT[1 << NUM_DIR] = { };

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
        offset_to_neighbor_even         [6] = 0;

        offset_to_eastern_boundary_even [0] = 0;
        offset_to_eastern_boundary_even [1] = 0;
        offset_to_eastern_boundary_even [2] = dim_x;
        offset_to_eastern_boundary_even [3] = dim_x;
        offset_to_eastern_boundary_even [4] = dim_x;
        offset_to_eastern_boundary_even [5] = 0;
        offset_to_eastern_boundary_even [6] = 0;

        offset_to_northern_boundary_even[0] = 0;
        offset_to_northern_boundary_even[1] = 0;
        offset_to_northern_boundary_even[2] = 0;
        offset_to_northern_boundary_even[3] = 0;
        offset_to_northern_boundary_even[4] = dim_x * dim_y;
        offset_to_northern_boundary_even[5] = dim_x * dim_y;
        offset_to_northern_boundary_even[6] = 0;

        offset_to_western_boundary_even [0] = -dim_x;
        offset_to_western_boundary_even [1] = 0;
        offset_to_western_boundary_even [2] = 0;
        offset_to_western_boundary_even [3] = 0;
        offset_to_western_boundary_even [4] = 0;
        offset_to_western_boundary_even [5] = 0;
        offset_to_western_boundary_even [6] = 0;

        offset_to_southern_boundary_even[0] = 0;
        offset_to_southern_boundary_even[1] = -dim_x * dim_y;
        offset_to_southern_boundary_even[2] = -dim_x * dim_y + 1;
        offset_to_southern_boundary_even[3] = 0;
        offset_to_southern_boundary_even[4] = 0;
        offset_to_southern_boundary_even[5] = 0;
        offset_to_southern_boundary_even[6] = 0;

        // Cell located in a row with odd index value
        offset_to_neighbor_odd          [0] = 1;
        offset_to_neighbor_odd          [1] = dim_x + 1;
        offset_to_neighbor_odd          [2] = dim_x;
        offset_to_neighbor_odd          [3] = -1;
        offset_to_neighbor_odd          [4] = -dim_x;
        offset_to_neighbor_odd          [5] = -dim_x + 1;
        offset_to_neighbor_odd          [6] = 0;

        offset_to_eastern_boundary_odd  [0] = 0;
        offset_to_eastern_boundary_odd  [1] = 0;
        offset_to_eastern_boundary_odd  [2] = 0;
        offset_to_eastern_boundary_odd  [3] = dim_x;
        offset_to_eastern_boundary_odd  [4] = 0;
        offset_to_eastern_boundary_odd  [5] = 0;
        offset_to_eastern_boundary_odd  [6] = 0;

        offset_to_northern_boundary_odd [0] = 0;
        offset_to_northern_boundary_odd [1] = 0;
        offset_to_northern_boundary_odd [2] = 0;
        offset_to_northern_boundary_odd [3] = 0;
        offset_to_northern_boundary_odd [4] = dim_x * dim_y;
        offset_to_northern_boundary_odd [5] = dim_x * dim_y;
        offset_to_northern_boundary_odd [6] = 0;

        offset_to_western_boundary_odd  [0] = -dim_x;
        offset_to_western_boundary_odd  [1] = -dim_x;
        offset_to_western_boundary_odd  [2] = 0;
        offset_to_western_boundary_odd  [3] = 0;
        offset_to_western_boundary_odd  [4] = 0;
        offset_to_western_boundary_odd  [5] = -dim_x;
        offset_to_western_boundary_odd  [6] = 0;

        offset_to_southern_boundary_odd [0] = 0;
        offset_to_southern_boundary_odd [1] = -dim_x * dim_y;
        offset_to_southern_boundary_odd [2] = -dim_x * dim_y;
        offset_to_southern_boundary_odd [3] = 0;
        offset_to_southern_boundary_odd [4] = 0;
        offset_to_southern_boundary_odd [5] = 0;
        offset_to_southern_boundary_odd [6] = 0;
    }

    static inline void collide(unsigned char* node_state_in, unsigned char* node_state_out, const bool p)
    {
        unsigned char a = node_state_in[1];
        unsigned char b = node_state_in[2];
        unsigned char c = node_state_in[3];
        unsigned char d = node_state_in[4];
        unsigned char e = node_state_in[5];
        unsigned char f = node_state_in[0];
        unsigned char r = node_state_in[6];

        unsigned char db1 = (a&d&~(b|c|e|f));
        unsigned char db2 = (b&e&~(a|c|d|f));
        unsigned char db3 = (c&f&~(a|b|d|e));

        unsigned char triple = (a^b)&(b^c)&(c^d)&(d^e)&(e^f);

        unsigned char ra = (r&a&~(b|c|d|e|f));
        unsigned char rb = (r&b&~(a|c|d|e|f));
        unsigned char rc = (r&c&~(a|b|d|e|f));
        unsigned char rd = (r&d&~(a|b|c|e|f));
        unsigned char re = (r&e&~(a|b|c|d|f));
        unsigned char rf = (r&f&~(a|b|c|d|e));

        unsigned char ra2 = (f&b&~(r|a|c|d|e));
        unsigned char rb2 = (a&c&~(r|b|d|e|f));
        unsigned char rc2 = (b&d&~(r|a|c|e|f));
        unsigned char rd2 = (c&e&~(r|a|b|d|f));
        unsigned char re2 = (d&f&~(r|a|b|c|e));
        unsigned char rf2 = (e&a&~(r|b|c|d|f));

        unsigned char adbe = db1&db2&(~(c|f));
        unsigned char adcf = db1&db3&(~(b|e));
        unsigned char becf = db2&db3&(~(a|d));

        unsigned char chad = (p&adbe)|((~p)&adcf)|becf;
        unsigned char chbe = (p&becf)|((~p)&adbe)|adcf;
        unsigned char chcf = (p&adcf)|((~p)&becf)|adbe;

        unsigned char adb = db1&b&(~(c|e|f));
        unsigned char adc = db1&c&(~(b|e|f));
        unsigned char ade = db1&e&(~(b|c|f));
        unsigned char adf = db1&f&(~(b|c|e));
        unsigned char bea = db2&a&(~(c|d|f));
        unsigned char bec = db2&c&(~(a|d|f));
        unsigned char bed = db2&d&(~(a|c|f));
        unsigned char bef = db2&f&(~(a|c|d));
        unsigned char cfa = db3&a&(~(b|d|e));
        unsigned char cfb = db3&b&(~(a|d|e));
        unsigned char cfd = db3&d&(~(a|b|e));
        unsigned char cfe = db3&e&(~(a|b|d));

        unsigned char spchad = (adb|ade|adc|adf)|(bec|bef)|(cfb|cfe);
        unsigned char spchbe = (bea|bec|bed|bef)|(adc|adf)|(cfa|cfd);
        unsigned char spchcf = (cfa|cfb|cfd|cfe)|(adb|ade)|(bea|bed);

        unsigned char cha = ((triple|db1|(p&db2)|((~p)&db3)|ra|rb|rf|ra2|rb2|rf2|spchad|chad)&(~0));
        unsigned char chd = ((triple|db1|(p&db2)|((~p)&db3)|rd|rc|re|rd2|rc2|re2|spchad|chad)&(~0));
        unsigned char chb = ((triple|db2|(p&db3)|((~p)&db1)|rb|ra|rc|rb2|ra2|rc2|spchbe|chbe)&(~0));
        unsigned char che = ((triple|db2|(p&db3)|((~p)&db1)|re|rd|rf|re2|rd2|rf2|spchbe|chbe)&(~0));
        unsigned char chc = ((triple|db3|(p&db1)|((~p)&db2)|rc|rb|rd|rc2|rb2|rd2|spchcf|chcf)&(~0));
        unsigned char chf = ((triple|db3|(p&db1)|((~p)&db2)|rf|ra|re|rf2|ra2|re2|spchcf|chcf)&(~0));
        unsigned char chr = ((ra|rb|rc|rd|re|rf|ra2|rb2|rc2|rd2|re2|rf2)&(~0));

        node_state_out[1] = node_state_in[1]^cha;
        node_state_out[2] = node_state_in[2]^chb;
        node_state_out[3] = node_state_in[3]^chc;
        node_state_out[4] = node_state_in[4]^chd;
        node_state_out[5] = node_state_in[5]^che;
        node_state_out[0] = node_state_in[0]^chf;
        node_state_out[6] = node_state_in[6]^chr;
    }

    static inline void bounce_back(unsigned char* node_state_in, unsigned char* node_state_out)
    {
#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

            // Exchange the states of the nodes with the the states of the inverse directions
            node_state_out[dir] = node_state_in[INV_DIR[dir]];
        }
    }

    static inline void bounce_forward_x(unsigned char* node_state_in, unsigned char* node_state_out)
    {
#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

                // Exchange the states of the nodes with the the states of the mirrored
                // directions along the x axis
                node_state_out[dir] = node_state_in[MIR_DIR_X[dir]];
        }
    }

    static inline void bounce_forward_y(unsigned char* node_state_in, unsigned char* node_state_out)
    {
#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir) {

                // Exchange the states of the nodes with the the states of
                // the mirrored directions along the y axis
                node_state_out[dir] = node_state_in[MIR_DIR_Y[dir]];
        }
    }

}; // struct ModelDescriptor<Model::FHP_III>

} // namespace lgca

#endif /* LGCA_MODELS_H_ */
