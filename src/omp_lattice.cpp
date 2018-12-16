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

#include "lgca_common.h"

#include "omp_lattice.h"

#include <omp.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace lgca {

// Definitions of static members
constexpr char          ModelDescriptor<Model::HPP>::INV_DIR[];
constexpr char          ModelDescriptor<Model::HPP>::MIR_DIR_X[];
constexpr char          ModelDescriptor<Model::HPP>::MIR_DIR_Y[];
constexpr Real          ModelDescriptor<Model::HPP>::LATTICE_VEC_X[];
constexpr Real          ModelDescriptor<Model::HPP>::LATTICE_VEC_Y[];

constexpr char          ModelDescriptor<Model::FHP_I>::INV_DIR[];
constexpr char          ModelDescriptor<Model::FHP_I>::MIR_DIR_X[];
constexpr char          ModelDescriptor<Model::FHP_I>::MIR_DIR_Y[];
constexpr Real          ModelDescriptor<Model::FHP_I>::LATTICE_VEC_X[];
constexpr Real          ModelDescriptor<Model::FHP_I>::LATTICE_VEC_Y[];

constexpr char          ModelDescriptor<Model::FHP_II>::INV_DIR[];
constexpr char          ModelDescriptor<Model::FHP_II>::MIR_DIR_X[];
constexpr char          ModelDescriptor<Model::FHP_II>::MIR_DIR_Y[];
constexpr Real          ModelDescriptor<Model::FHP_II>::LATTICE_VEC_X[];
constexpr Real          ModelDescriptor<Model::FHP_II>::LATTICE_VEC_Y[];

constexpr char          ModelDescriptor<Model::FHP_III>::INV_DIR[];
constexpr char          ModelDescriptor<Model::FHP_III>::MIR_DIR_X[];
constexpr char          ModelDescriptor<Model::FHP_III>::MIR_DIR_Y[];
constexpr Real          ModelDescriptor<Model::FHP_III>::LATTICE_VEC_X[];
constexpr Real          ModelDescriptor<Model::FHP_III>::LATTICE_VEC_Y[];


// Creates a CUDA parallelized lattice gas cellular automaton object
// of the specified properties.
template<Model model_>
OMP_Lattice<model_>::OMP_Lattice(const string test_case,
                                 const Real Re, const Real Ma_s,
                                 const int coarse_graining_radius)
               : Lattice<model_>(test_case, Re, Ma_s, coarse_graining_radius) {

    // Allocate the memory for the arrays on the host (CPU)
    allocate_memory();

    // Generate random bits for collision
    this->m_rnd_cpu.fill_random();

    // Set the model-based values according to the number of lattice directions
    m_model = new ModelDesc(this->m_dim_x, this->m_dim_y);
}

// Deletes the openMP parallelized lattice gas cellular automaton object.
template<Model model_>
OMP_Lattice<model_>::~OMP_Lattice() {

    delete m_model;

    this->free_memory();
}

// Performs the collision and propagation step on the lattice gas automaton.
template<Model model_>
void OMP_Lattice<model_>::collide_and_propagate() {

    Bitset::Block* __restrict__ read  = this->m_node_state_cpu.ptr();
    Bitset::Block* __restrict__ write = this->m_node_state_tmp_cpu.ptr();

    // Loop over bunches of cells
    const size_t num_cell_blocks = ((this->m_num_cells - 1) / Bitset::BITS_PER_BLOCK) + 1;
    const size_t cell_block_dim_x = this->m_dim_x / Bitset::BITS_PER_BLOCK;
    // TODO Parallelization leads to incorrect results!
//    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_cell_blocks), [&](const tbb::blocked_range<size_t>& r) {
//    for (size_t cell_block = r.begin(); cell_block != r.end(); ++cell_block)
    for (size_t cell_block = 0; cell_block < num_cell_blocks; ++cell_block)
    {
        const size_t cell_block_pos_x = cell_block % cell_block_dim_x;
        const size_t cell_block_pos_y = cell_block / cell_block_dim_x;

        if (cell_block_pos_x == 0 || cell_block_pos_x == cell_block_dim_x-1 ||
            cell_block_pos_y == 0 || cell_block_pos_y == this->m_dim_y-1) continue;

        Bitset::Block is_fluid    = this->m_cell_type_cpu(cell_block);
        Bitset::Block random_bits = this->m_rnd_cpu(cell_block);

        Bitset::Block inputs[this->NUM_DIR];
#pragma unroll
        for (int dir = 0; dir < this->NUM_DIR; ++dir)
            inputs[dir] = read[cell_block * this->NUM_DIR + dir];

        // By default, transport
        Bitset::Block outputs[this->NUM_DIR];
#pragma unroll
        for (int dir = 0; dir < this->NUM_DIR; ++dir)
            outputs[dir] = inputs[dir];

        // Model head on collisions
        const Bitset::Block collision0 =  inputs[0] & ~inputs[1] & ~inputs[2] &  inputs[3] & ~inputs[4] & ~inputs[5] & is_fluid;
        const Bitset::Block collision1 = ~inputs[0] &  inputs[1] & ~inputs[2] & ~inputs[3] &  inputs[4] & ~inputs[5] & is_fluid;
        const Bitset::Block collision2 = ~inputs[0] & ~inputs[1] &  inputs[2] & ~inputs[3] & ~inputs[4] &  inputs[5] & is_fluid;

        outputs[0] &= ~collision0; outputs[3] &= ~collision0;
        outputs[1] |=  collision0 &  random_bits; outputs[4] |= collision0 &  random_bits;
        outputs[2] |=  collision0 & ~random_bits; outputs[5] |= collision0 & ~random_bits;

        outputs[1] &= ~collision1; outputs[4] &= ~collision1;
        outputs[2] |= collision1 &  random_bits; outputs[5] |= collision1 &  random_bits;
        outputs[3] |= collision1 & ~random_bits; outputs[0] |= collision1 & ~random_bits;

        outputs[2] &= ~collision2; outputs[5] &= ~collision2;
        outputs[3] |= collision2 &  random_bits; outputs[0] |= collision2 &  random_bits;
        outputs[4] |= collision2 & ~random_bits; outputs[1] |= collision2 & ~random_bits;

        // Model three way collisions
        const Bitset::Block collision3 =  inputs[0] & ~inputs[1] &  inputs[2] & ~inputs[3] &  inputs[4] & ~inputs[5] & is_fluid;
        const Bitset::Block collision4 = ~inputs[0] &  inputs[1] & ~inputs[2] &  inputs[3] & ~inputs[4] &  inputs[5] & is_fluid;

        outputs[0] &= ~collision3; outputs[2] &= ~collision3; outputs[4] &= ~collision3;
        outputs[1] |=  collision3; outputs[3] |=  collision3; outputs[5] |=  collision3;

        outputs[1] &= ~collision4; outputs[3] &= ~collision4; outputs[5] &= ~collision4;
        outputs[0] |=  collision4; outputs[2] |=  collision4; outputs[4] |=  collision4;

        // Model no slip boundary conditions (bounce back)
        outputs[0] ^= (outputs[0] ^ inputs[3]) & ~is_fluid;
        outputs[1] ^= (outputs[1] ^ inputs[4]) & ~is_fluid;
        outputs[2] ^= (outputs[2] ^ inputs[5]) & ~is_fluid;
        outputs[3] ^= (outputs[3] ^ inputs[0]) & ~is_fluid;
        outputs[4] ^= (outputs[4] ^ inputs[1]) & ~is_fluid;
        outputs[5] ^= (outputs[5] ^ inputs[2]) & ~is_fluid;

        // TODO Adaptive w.r.t. bitset data type
        Bitset::Block odd_col = 6148914691236517205; // binary 1010101010101010...10

        // Push outputs (all columns)
        write[(cell_block_pos_x + (cell_block_pos_y+1) * cell_block_dim_x) * this->NUM_DIR + 0]  = outputs[0];
        write[(cell_block_pos_x + (cell_block_pos_y-1) * cell_block_dim_x) * this->NUM_DIR + 3]  = outputs[3];

        // Push outputs (even columns)
        Bitset::Block* x0y01 = &write[(cell_block_pos_x   + (cell_block_pos_y  ) * cell_block_dim_x) * this->NUM_DIR + 1];
        Bitset::Block* xmy01 = &write[(cell_block_pos_x-1 + (cell_block_pos_y  ) * cell_block_dim_x) * this->NUM_DIR + 1];
        Bitset::Block* x0ym2 = &write[(cell_block_pos_x   + (cell_block_pos_y-1) * cell_block_dim_x) * this->NUM_DIR + 2];
        Bitset::Block* xmym2 = &write[(cell_block_pos_x-1 + (cell_block_pos_y-1) * cell_block_dim_x) * this->NUM_DIR + 2];
        Bitset::Block* x0ym4 = &write[(cell_block_pos_x   + (cell_block_pos_y-1) * cell_block_dim_x) * this->NUM_DIR + 4];
        Bitset::Block* xpym4 = &write[(cell_block_pos_x+1 + (cell_block_pos_y-1) * cell_block_dim_x) * this->NUM_DIR + 4];
        Bitset::Block* x0y05 = &write[(cell_block_pos_x   +  cell_block_pos_y    * cell_block_dim_x) * this->NUM_DIR + 5];
        Bitset::Block* xpy05 = &write[(cell_block_pos_x+1 +  cell_block_pos_y    * cell_block_dim_x) * this->NUM_DIR + 5];

        *x0y01 ^= ((*x0y01) ^ (outputs[1] >> 1))                                     & ~odd_col;
        *xmy01 ^= ((*xmy01) ^ ((*xmy01) | (outputs[1] << Bitset::BITS_PER_BLOCK-1))) & ~odd_col;
        *x0ym2 ^= ((*x0ym2) ^ ((*x0ym2) | (outputs[2] >> 1)))                        & ~odd_col;
        *xmym2 ^= ((*xmym2) ^ ((*xmym2) | (outputs[2] << Bitset::BITS_PER_BLOCK-1))) & ~odd_col;
        *x0ym4 ^= ((*x0ym4) ^ ((*x0ym4) | (outputs[4] << 1)))                        & ~odd_col;
        *xpym4 ^= ((*xpym4) ^ ((*xpym4) | (outputs[4] >> Bitset::BITS_PER_BLOCK-1))) & ~odd_col;
        *x0y05 ^= ((*x0y05) ^ ((*x0y05) | (outputs[5] << 1)))                        & ~odd_col;
        *xpy05 ^= ((*xpy05) ^ (outputs[5] >> Bitset::BITS_PER_BLOCK-1))              & ~odd_col;

        // Push outputs (odd columns)
        Bitset::Block* x0yp1 = &write[(cell_block_pos_x   + (cell_block_pos_y+1) * cell_block_dim_x) * this->NUM_DIR + 1];
        Bitset::Block* xmyp1 = &write[(cell_block_pos_x-1 + (cell_block_pos_y+1) * cell_block_dim_x) * this->NUM_DIR + 1];
        Bitset::Block* x0y02 = &write[(cell_block_pos_x   +  cell_block_pos_y    * cell_block_dim_x) * this->NUM_DIR + 2];
        Bitset::Block* xmy02 = &write[(cell_block_pos_x-1 +  cell_block_pos_y    * cell_block_dim_x) * this->NUM_DIR + 2];
        Bitset::Block* x0y04 = &write[(cell_block_pos_x   +  cell_block_pos_y    * cell_block_dim_x) * this->NUM_DIR + 4];
        Bitset::Block* xpy04 = &write[(cell_block_pos_x+1 +  cell_block_pos_y    * cell_block_dim_x) * this->NUM_DIR + 4];
        Bitset::Block* x0yp5 = &write[(cell_block_pos_x   + (cell_block_pos_y+1) * cell_block_dim_x) * this->NUM_DIR + 5];
        Bitset::Block* xpyp5 = &write[(cell_block_pos_x+1 + (cell_block_pos_y+1) * cell_block_dim_x) * this->NUM_DIR + 5];

        *x0yp1 ^= ((*x0yp1) ^ (outputs[1] >> 1))                                     & odd_col;
        *xmyp1 ^= ((*xmyp1) ^ ((*xmyp1) | (outputs[1] << Bitset::BITS_PER_BLOCK-1))) & odd_col;
        *x0y02 ^= ((*x0y02) ^ ((*x0y02) | (outputs[2] >> 1)))                        & odd_col;
        *xmy02 ^= ((*xmy02) ^ ((*xmy02) | (outputs[2] << Bitset::BITS_PER_BLOCK-1))) & odd_col;
        *x0y04 ^= ((*x0y04) ^ ((*x0y04) | (outputs[4] << 1)))                        & odd_col;
        *xpy04 ^= ((*xpy04) ^ ((*xpy04) | (outputs[4] >> Bitset::BITS_PER_BLOCK-1))) & odd_col;
        *x0yp5 ^= ((*x0yp5) ^ ((*x0yp5) | (outputs[5] << 1)))                        & odd_col;
        *xpyp5 ^= ((*xpyp5) ^ (outputs[5] >> Bitset::BITS_PER_BLOCK-1))              & odd_col;

    }//}); // for cell block

    // Enforce periodic boundary conditions
    //
    // Top and bottom rows
    for (size_t xx = 0; xx < cell_block_dim_x; xx++) {
        for (int dir = 0; dir < this->NUM_DIR; dir++) {
            write[(xx + (this->m_dim_y-2) * cell_block_dim_x) * this->NUM_DIR + dir] ^=
            write[ xx                                         * this->NUM_DIR + dir];
        }
        for (int dir = 0; dir < this->NUM_DIR; dir++) {
            write[(xx +                     cell_block_dim_x) * this->NUM_DIR + dir] ^=
            write[(xx + (this->m_dim_y-1) * cell_block_dim_x) * this->NUM_DIR + dir];
        }
    }

    // Left and right sides
    for (size_t yy = 0; yy < this->m_dim_y; yy++) {
        for (int dir = 0; dir < this->NUM_DIR; dir++) {
            write[(cell_block_dim_x-2 + yy * cell_block_dim_x) * this->NUM_DIR + dir] ^=
            write[(                     yy * cell_block_dim_x) * this->NUM_DIR + dir];
        }
        for (int dir = 0; dir < this->NUM_DIR; dir++) {
            write[(                 1 + yy * cell_block_dim_x) * this->NUM_DIR + dir] ^=
            write[(cell_block_dim_x-1 + yy * cell_block_dim_x) * this->NUM_DIR + dir];
        }
    }

    this->m_node_state_cpu.reset();

    // Update the node states
    auto node_state_cpu_tmp = this->m_node_state_cpu.ptr();
    this->m_node_state_cpu  = m_node_state_tmp_cpu.ptr();
    m_node_state_tmp_cpu = node_state_cpu_tmp;
}

// Applies a body force in the specified direction (x or y) and with the
// specified intensity to the particles. E.g., if the intensity is equal 100,
// every 100th particle changes it's direction, if feasible.
template<Model model_>
void OMP_Lattice<model_>::apply_body_force(const int forcing) {

    // Set a maximum number of iterations to find particles which can be reverted
    const size_t it_max = 2 * this->m_num_cells; size_t it = 0;

    // Number of particles which have been reverted
    unsigned int reverted_particles = 0;

    // Loop over all cells
    do
    {
        it++;

        // Choose a random cell
        size_t global_cell = rand() % this->m_num_cells;

        // Get the type of the cell, i.e. fluid or solid.
        // Note that body forces are applied to fluid cells only.
        bool cell_type = bool(this->m_cell_type_cpu[global_cell]);

        // Check weather the cell working on is a fluid cell
        if (cell_type == bool(CellType::FLUID)) {

            // Get the whole cell block containing the cell under consideration
            Bitset node_state(Bitset::BITS_PER_BLOCK * this->NUM_DIR);
            const size_t cell_block = global_cell / Bitset::BITS_PER_BLOCK;
            const size_t local_cell = global_cell % Bitset::BITS_PER_BLOCK;
#pragma unroll
            for (int dir = 0; dir < this->NUM_DIR; ++dir)
                node_state(dir) = this->m_node_state_cpu(cell_block * this->NUM_DIR + dir);

            if (model_ == Model::HPP) {

//                if (this->m_bf_dir == 'x' && (node_state[0] == 0) && (node_state[2] == 1)) {

//					node_state_tmp[0] = 1;
//					node_state_tmp[2] = 0;

//					reverted_particles++;

//                } else if (this->m_bf_dir == 'y' && (node_state[1] == 1) && (node_state[3] == 0)) {

//					node_state_tmp[1] = 0;
//					node_state_tmp[3] = 1;

//					reverted_particles++;
//				}
        	}

            else if (model_ == Model::FHP_I || model_ == Model::FHP_II || model_ == Model::FHP_III) {

                if (this->m_bf_dir == 'x' && (node_state[5 * Bitset::BITS_PER_BLOCK + local_cell] == bool(0))
                                          && (node_state[1 * Bitset::BITS_PER_BLOCK + local_cell] == bool(1))) {

                    node_state[5 * Bitset::BITS_PER_BLOCK + local_cell] = bool(1);
                    node_state[1 * Bitset::BITS_PER_BLOCK + local_cell] = bool(0);

					reverted_particles++;

                } else if (this->m_bf_dir == 'y') {

//                    if ((node_state[1] == 1) && (node_state[5] == 0)) {

//                        node_state_tmp[1] = 0;
//                        node_state_tmp[5] = 1;

//                        reverted_particles++;
//                    }

//                    if ((node_state[2] == 1) && (node_state[4] == 0)) {

//                        node_state_tmp[2] = 0;
//                        node_state_tmp[4] = 1;

//                        reverted_particles++;
//                    }
                }
        	}

            // Write the new node states back to the data array
#pragma unroll
            for (int dir = 0; dir < this->NUM_DIR; ++dir)
                this->m_node_state_cpu(cell_block * this->NUM_DIR + dir) = node_state(dir);

        } /* IF cell_type */

    } while ((reverted_particles < forcing) && (it < it_max));
}

// Computes quantities of interest as a post-processing procedure
template<Model model_>
void OMP_Lattice<model_>::post_process() {

    // Computes cell quantities of interest as a post-processing procedure
	cell_post_process();

    // Computes coarse grained quantities of interest as a post-processing procedure
	mean_post_process();
}

// Computes cell quantities of interest as a post-processing procedure
template<Model model_>
void OMP_Lattice<model_>::cell_post_process()
{
    const Bitset::Block* __restrict__ read = this->m_node_state_out_cpu.ptr();

    // Loop over bunches of cells
    const size_t num_cell_blocks = ((this->m_num_cells - 1) / Bitset::BITS_PER_BLOCK) + 1;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_cell_blocks), [&](const tbb::blocked_range<size_t>& r) {
    for (size_t cell_block = r.begin(); cell_block != r.end(); ++cell_block)
    {
        Bitset inputs(Bitset::BITS_PER_BLOCK * this->NUM_DIR);
#pragma unroll
        for (int dir = 0; dir < this->NUM_DIR; ++dir)
            inputs(dir) = read[cell_block * this->NUM_DIR + dir];

        // Loop over cells in cell block
#pragma unroll
        for (int local_cell = 0; local_cell < Bitset::BITS_PER_BLOCK; ++local_cell) {

            // Initialize the cell quantities to be computed
            unsigned char cell_density    = 0;
            Real          cell_momentum_x = 0.0;
            Real          cell_momentum_y = 0.0;
            Real          cell_momentum_z = 0.0;

            // Loop over nodes within the current cell
#pragma unroll
            for (int dir = 0; dir < this->NUM_DIR; ++dir) {

                bool node_state = bool(inputs[local_cell + dir * Bitset::BITS_PER_BLOCK]);

                // Sum up the node states
                cell_density += node_state;

                // Sum up the node states multiplied by the lattice vector component for the current
                // direction
                cell_momentum_x += node_state * ModelDesc::LATTICE_VEC_X[dir];
                cell_momentum_y += node_state * ModelDesc::LATTICE_VEC_Y[dir];
            }

            // Write the computed cell quantities to the related data arrays
            const size_t global_cell = cell_block * Bitset::BITS_PER_BLOCK + local_cell;
            this->m_cell_density_cpu [global_cell] = (Real) cell_density;

            Real cell_vel_x = cell_momentum_x / Real(cell_density > 0 ? cell_density : 1);
            Real cell_vel_y = cell_momentum_y / Real(cell_density > 0 ? cell_density : 1);
            Real cell_vel_z = cell_momentum_z / Real(cell_density > 0 ? cell_density : 1);

            this->m_cell_velocity_cpu[global_cell * 3 + 0] = cell_vel_x;
            this->m_cell_velocity_cpu[global_cell * 3 + 1] = cell_vel_y;
            this->m_cell_velocity_cpu[global_cell * 3 + 2] = cell_vel_z;

        } // for local cell

    }}); // for cell block
}

// Computes coarse grained quantities of interest as a post-processing procedure
template<Model model_>
void OMP_Lattice<model_>::mean_post_process()
{
    const int r = this->m_coarse_graining_radius;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, this->m_num_coarse_cells), [&](const tbb::blocked_range<size_t>& range) {
    for (size_t coarse_cell = range.begin(); coarse_cell != range.end(); ++coarse_cell)
    {
        // Get cell in the bottom left corner of the coarse cell
        const size_t cell = (coarse_cell % this->m_coarse_dim_x) * (2 * r)
                          + (coarse_cell / this->m_coarse_dim_x) * (2 * r) * this->m_dim_x;

        // Calculate the position of the cell in x direction
        int pos_x = cell % this->m_dim_x;

        // Initialize the coarse grained quantities to be computed
        Real mean_density    = 0.0;
        Real mean_momentum_x = 0.0;
        Real mean_momentum_y = 0.0;
        Real mean_momentum_z = 0.0;

        // Initialize the number of actual existing coarse graining neighbor cells
        int n_exist_neighbors = 0;

        // The thread working on the cell has to know the cell quantities of the coarse graining
        // neighbor cells, therefore looping over all neighbor cells and look it up
#pragma unroll
        for (int y = 0; y <= 2 * r; ++y) {

            for (int x = 0; x <= 2 * r; ++x) {

                // Get the index of the coarse graining neighbor cell
                size_t neighbor_idx = cell + y * this->m_dim_x + x;

                // Get the position of the coarse graining neighbor cell in x direction
                int pos_x_neighbor = neighbor_idx % this->m_dim_x;

                // Check weather the coarse graining neighbor cell is valid
                if ((neighbor_idx >= 0) &&
                    (neighbor_idx < this->m_num_cells) &&
                    (abs(pos_x_neighbor - pos_x) <= r)) {

                    // Increase the number of existing coarse graining neighbor cells
                    n_exist_neighbors++;

                    mean_density    += this->m_cell_density_cpu [neighbor_idx        ];
                    mean_momentum_x += this->m_cell_velocity_cpu[neighbor_idx * 3 + 0];
                    mean_momentum_y += this->m_cell_velocity_cpu[neighbor_idx * 3 + 1];
                }
            }
        }

        // Write the computed coarse grained quantities to the related data arrays
        this->m_mean_density_cpu [coarse_cell        ] = mean_density    / ((Real) n_exist_neighbors);
        this->m_mean_velocity_cpu[coarse_cell * 3 + 0] = mean_momentum_x / ((Real) n_exist_neighbors);
        this->m_mean_velocity_cpu[coarse_cell * 3 + 1] = mean_momentum_y / ((Real) n_exist_neighbors);
        this->m_mean_velocity_cpu[coarse_cell * 3 + 2] = mean_momentum_z / ((Real) n_exist_neighbors);

    }}); // for coarse_cell
}

// Allocates the memory for the arrays on the host (CPU)
template<Model model_>
void OMP_Lattice<model_>::allocate_memory()
{
    // Allocate host memory
    this->m_cell_density_cpu   = (Real*)malloc(    this->m_num_cells        * sizeof(Real));
    this->m_mean_density_cpu   = (Real*)malloc(    this->m_num_coarse_cells * sizeof(Real));
    this->m_cell_velocity_cpu  = (Real*)malloc(3 * this->m_num_cells        * sizeof(Real));
    this->m_mean_velocity_cpu  = (Real*)malloc(3 * this->m_num_coarse_cells * sizeof(Real));

    this->m_node_state_cpu.resize    (this->m_num_nodes);
    this->m_node_state_tmp_cpu.resize(this->m_num_nodes);
    this->m_node_state_out_cpu.resize(this->m_num_nodes);
    this->m_rnd_cpu.resize           (this->m_num_cells);
    this->m_cell_type_cpu.resize     (this->m_num_cells);
}

// Frees the memory for the arrays on the host (CPU)
template<Model model_>
void OMP_Lattice<model_>::free_memory()
{
    // Free CPU memory
    free(this->m_cell_density_cpu);
    free(this->m_mean_density_cpu);
    free(this->m_cell_velocity_cpu);
    free(this->m_mean_velocity_cpu);

    this->m_cell_density_cpu    = NULL;
    this->m_mean_density_cpu    = NULL;
    this->m_cell_velocity_cpu   = NULL;
    this->m_mean_velocity_cpu   = NULL;
}

// Computes the mean velocity of the lattice.
template<Model model_>
std::vector<Real> OMP_Lattice<model_>::get_mean_velocity() {

    std::vector<Real> mean_velocity(this->SPATIAL_DIM, 0.0);

    Real sum_x_vel = 0.0;
    Real sum_y_vel = 0.0;

    size_t counter = 0;

    // Sum up all (fluid) cell x and y velocity components.
#pragma omp parallel for reduction(+: sum_x_vel, sum_y_vel)
    for (size_t n = 0; n < this->m_num_cells; ++n) {

        if (this->m_cell_type_cpu[n] == bool(CellType::FLUID)) {

        	counter++;

            Real cell_density = this->m_cell_density_cpu[n];

			if (cell_density > 1.0e-06) {

                sum_x_vel += this->m_cell_velocity_cpu[n * 3 + 0] / cell_density;
                sum_y_vel += this->m_cell_velocity_cpu[n * 3 + 1] / cell_density;
			}

#ifdef DEBUG

			else if (fabs(cell_density) < 1.0e-06) {

				// Do nothing.

			} else if (cell_density < -1.0e-06) {

				printf("ERROR in get_mean_velocity(): "
					   "Negative cell density detected.");
				abort();
			}

#endif

        }
    }

    // Divide the summed up x and y components by the total number of fluid cells.
    mean_velocity[0] = sum_x_vel / (Real) counter;
    mean_velocity[1] = sum_y_vel / (Real) counter;

    return mean_velocity;
}

// Explicit instantiations
template class OMP_Lattice<Model::HPP>;
template class OMP_Lattice<Model::FHP_I>;
template class OMP_Lattice<Model::FHP_II>;
template class OMP_Lattice<Model::FHP_III>;

} // namespace lgca
