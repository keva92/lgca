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
#include "cuda_utils.cuh"

#include <omp.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace lgca {

// Definitions of static members
constexpr char ModelDescriptor<Model::HPP>::INV_DIR[];
constexpr char ModelDescriptor<Model::HPP>::MIR_DIR_X[];
constexpr char ModelDescriptor<Model::HPP>::MIR_DIR_Y[];
constexpr Real ModelDescriptor<Model::HPP>::LATTICE_VEC_X[];
constexpr Real ModelDescriptor<Model::HPP>::LATTICE_VEC_Y[];

constexpr char ModelDescriptor<Model::FHP>::INV_DIR[];
constexpr char ModelDescriptor<Model::FHP>::MIR_DIR_X[];
constexpr char ModelDescriptor<Model::FHP>::MIR_DIR_Y[];
constexpr Real ModelDescriptor<Model::FHP>::LATTICE_VEC_X[];
constexpr Real ModelDescriptor<Model::FHP>::LATTICE_VEC_Y[];


// Creates a CUDA parallelized lattice gas cellular automaton object
// of the specified properties.
template<Model model_>
OMP_Lattice<model_>::OMP_Lattice(const string test_case,
                                 const Real Re, const Real Ma_s,
                                 const int coarse_graining_radius)
               : Lattice<model_>(test_case, Re, Ma_s, coarse_graining_radius) {

    // Allocate the memory for the arrays on the host (CPU)
    allocate_memory();

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
void OMP_Lattice<model_>::collide_and_propagate(unsigned int step) {

    // TODO Set the seed for the random number generation.
    int seed = time(NULL);

#ifndef NDEBUG
            // Check weather the domain dimensions are valid for the FHP model.
            if (this->m_dim_y % 2 != 0 && model_ == Model::FHP) {

                printf("ERROR in collide_and_propagate_kernel(): "
                       "Invalid domain dimension in y direction.\n");
                abort();
            }
#endif

    // Loop over bunches of cells
    const int num_blocks = ((this->m_num_cells - 1) / Bitset::BITS_PER_BLOCK) + 1;
    tbb::parallel_for(tbb::blocked_range<int>(0, num_blocks), [&](const tbb::blocked_range<int>& r) {
    for (int block = r.begin(); block != r.end(); ++block)
    {        
        for (int cell = block * Bitset::BITS_PER_BLOCK; cell < (block+1) * Bitset::BITS_PER_BLOCK; ++cell) {

            if (cell >= this->m_num_cells) break;

            // Calculate the position of the cell in x direction (column index).
            int pos_x = cell % this->m_dim_x;

            // Calculate the position of the cell in y direction (row index).
            int pos_y = cell / this->m_dim_x;

            // Get the type of the cell, i.e. fluid or solid.
            // This has to be taken into account during the collision step, where
            // cells behave different according to their type.
            char cell_type = this->m_cell_type_cpu[cell];

            // Check weather the cell is located on boundaries.
            bool on_eastern_boundary  = (cell + 1) % this->m_dim_x == 0;
            bool on_northern_boundary = cell >= (this->m_num_cells - this->m_dim_x);
            bool on_western_boundary  = cell % this->m_dim_x == 0;
            bool on_southern_boundary = cell < this->m_dim_x;

            // Define an array for the global indices of the nodes in the cell.
            int node_idx[this->NUM_DIR];

            // Define an array for the states of the nodes in the cell.
            char node_state[this->NUM_DIR];

            // Execute collision step.
            //
            // The thread working on the cell has to know about the states of the
            // nodes within the cell, therefore looping over all directions and
            // look it up.
#pragma unroll
            for (int dir = 0; dir < this->NUM_DIR; ++dir) {

                node_idx[dir] = dir + cell * this->NUM_DIR;
                node_state[dir] = bool(this->m_node_state_cpu[node_idx[dir]]);
            }

            // TODO Create a random boolean value for the collision step.
//            bool rand_bool = cu_random_bool(seed, cell);
//            bool rand_bool =       ((pos_x % 2) == (pos_y % 2))
//                             - 1 * ((pos_x % 2) == (pos_y % 2)) * (step % 2)
//                             + 1 * ((pos_x % 2) != (pos_y % 2)) * (step % 2);
//            bool rand_bool = true;
            bool rand_bool = false;

            // Create a temporary array to copy the node states into.
            char node_state_tmp[this->NUM_DIR];

            // Copy the actual states of the nodes to the temporary array.
#pragma unroll
            for (int dir = 0; dir < this->NUM_DIR; ++dir) {

                node_state_tmp[dir] = node_state[dir];
            }

            switch (cell_type) {

                // The cell working on is a fluid cell ("normal" collision).
                case 0:
                {
                    // Using the the HPP model.
                    if (model_ == Model::HPP) {

        //                    // Collision case 1.
        //                    if ((node_state[0] == 0) &&
        //                        (node_state[1] == 1) &&
        //                        (node_state[2] == 0) &&
        //                        (node_state[3] == 1)) {
        //
        //                        node_state_tmp[0] = 1;
        //                        node_state_tmp[1] = 0;
        //                        node_state_tmp[2] = 1;
        //                        node_state_tmp[3] = 0;
        //
        //                        break;
        //                    }
        //
        //                    // Collision case 2.
        //                    if ((node_state[0] == 1) &&
        //                        (node_state[1] == 0) &&
        //                        (node_state[2] == 1) &&
        //                        (node_state[3] == 0)) {
        //
        //                        node_state_tmp[0] = 0;
        //                        node_state_tmp[1] = 1;
        //                        node_state_tmp[2] = 0;
        //                        node_state_tmp[3] = 1;
        //
        //                        break;
        //                    }

                        node_state_tmp[0] = node_state[0]
                                - (node_state[0] * node_state[2] * (1 - node_state[1]) * (1 - node_state[3]))
                                + (node_state[1] * node_state[3] * (1 - node_state[0]) * (1 - node_state[2]));

                        node_state_tmp[1] = node_state[1]
                                - (node_state[1] * node_state[3] * (1 - node_state[0]) * (1 - node_state[2]))
                                + (node_state[0] * node_state[2] * (1 - node_state[1]) * (1 - node_state[3]));

                        node_state_tmp[2] = node_state[2]
                                - (node_state[0] * node_state[2] * (1 - node_state[1]) * (1 - node_state[3]))
                                + (node_state[1] * node_state[3] * (1 - node_state[0]) * (1 - node_state[2]));

                        node_state_tmp[3] = node_state[3]
                                - (node_state[1] * node_state[3] * (1 - node_state[0]) * (1 - node_state[2]))
                                + (node_state[0] * node_state[2] * (1 - node_state[1]) * (1 - node_state[3]));

                    // Collision cases of the FHP model.
                    } else if (model_ == Model::FHP) {

                        // Collision case a1.
                        if ((node_state[0] == 1) &&
                            (node_state[1] == 0) &&
                            (node_state[2] == 0) &&
                            (node_state[3] == 1) &&
                            (node_state[4] == 0) &&
                            (node_state[5] == 0)) {

                            node_state_tmp[0] = 0;
                            node_state_tmp[1] = rand_bool;
                            node_state_tmp[2] = 1 - node_state_tmp[1];
                            node_state_tmp[3] = 0;
                            node_state_tmp[4] = node_state_tmp[1];
                            node_state_tmp[5] = node_state_tmp[2];

                            break;
                        }

                        // Collision case a2.
                        if ((node_state[0] == 0) &&
                            (node_state[1] == 1) &&
                            (node_state[2] == 0) &&
                            (node_state[3] == 0) &&
                            (node_state[4] == 1) &&
                            (node_state[5] == 0)) {

                            node_state_tmp[0] = rand_bool;
                            node_state_tmp[1] = 0;
                            node_state_tmp[2] = 1 - node_state_tmp[0];
                            node_state_tmp[3] = node_state_tmp[0];
                            node_state_tmp[4] = 0;
                            node_state_tmp[5] = node_state_tmp[2];

                            break;
                        }

                        // Collision case a3.
                        if ((node_state[0] == 0) &&
                            (node_state[1] == 0) &&
                            (node_state[2] == 1) &&
                            (node_state[3] == 0) &&
                            (node_state[4] == 0) &&
                            (node_state[5] == 1)) {

                            node_state_tmp[0] = rand_bool;
                            node_state_tmp[1] = 1 - node_state_tmp[0];
                            node_state_tmp[2] = 0;
                            node_state_tmp[3] = node_state_tmp[0];
                            node_state_tmp[4] = node_state_tmp[1];
                            node_state_tmp[5] = 0;

                            break;
                        }

                        // Collision case b1.
                        if ((node_state[0] == 0) &&
                            (node_state[1] == 1) &&
                            (node_state[2] == 0) &&
                            (node_state[3] == 1) &&
                            (node_state[4] == 0) &&
                            (node_state[5] == 1)) {

                            node_state_tmp[0] = 1;
                            node_state_tmp[1] = 0;
                            node_state_tmp[2] = 1;
                            node_state_tmp[3] = 0;
                            node_state_tmp[4] = 1;
                            node_state_tmp[5] = 0;

                            break;
                        }

                        // Collision case b2.
                        if ((node_state[0] == 1) &&
                            (node_state[1] == 0) &&
                            (node_state[2] == 1) &&
                            (node_state[3] == 0) &&
                            (node_state[4] == 1) &&
                            (node_state[5] == 0)) {

                            node_state_tmp[0] = 0;
                            node_state_tmp[1] = 1;
                            node_state_tmp[2] = 0;
                            node_state_tmp[3] = 1;
                            node_state_tmp[4] = 0;
                            node_state_tmp[5] = 1;

                            break;
                        }

        //                    node_state_tmp[0] = node_state[0]
        //                            - (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5]))
        //                            + (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
        //                            + (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4])) * (1 - rand_bool)
        //                            - (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]))
        //                            + (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]));
        //
        //                    node_state_tmp[1] = node_state[1]
        //                            - (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5]))
        //                            + (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
        //                            + (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4])) * (1 - rand_bool)
        //                            - (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]))
        //                            + (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]));
        //
        //                    node_state_tmp[2] = node_state[2]
        //                            - (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4]))
        //                            + (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
        //                            + (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5])) * (1 - rand_bool)
        //                            - (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]))
        //                            + (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]));
        //
        //                    node_state_tmp[3] = node_state[3]
        //                            - (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5]))
        //                            + (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
        //                            + (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4])) * (1 - rand_bool)
        //                            - (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]))
        //                            + (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]));
        //
        //                    node_state_tmp[4] = node_state[4]
        //                            - (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5]))
        //                            + (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
        //                            + (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4])) * (1 - rand_bool)
        //                            - (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]))
        //                            + (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]));
        //
        //                    node_state_tmp[5] = node_state[5]
        //                            - (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4]))
        //                            + (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
        //                            + (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5])) * (1 - rand_bool)
        //                            - (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]))
        //                            + (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]));
                    }

                    break;
                }

                // The cell working on is a solid cell of bounce back type.
                case 1:
                {
                    // Loop over all directions.
#pragma unroll
                    for (int dir = 0; dir < this->NUM_DIR; ++dir) {

                        // Exchange the states of the nodes with the the states of
                        // the inverse directions.
                        node_state_tmp[dir] = node_state[ModelDesc::INV_DIR[dir]];
                    }

                    break;
                }

                // TODO: The cell working on is a solid cell of bounce forward type.
                case 2:
                {
                    // Loop over all directions.
#pragma unroll
                    for (int dir = 0; dir < this->NUM_DIR; ++dir) {

                        if (on_northern_boundary || on_southern_boundary) {

                            // Exchange the states of the nodes with the the states of
                            // the mirrored directions along the x axis.
                            node_state_tmp[dir] = node_state[ModelDesc::MIR_DIR_X[dir]];
                        }

                        if (on_eastern_boundary || on_western_boundary) {

                            // Exchange the states of the nodes with the the states of
                            // the mirrored directions along the y axis.
                            node_state_tmp[dir] = node_state[ModelDesc::MIR_DIR_Y[dir]];
                        }
                    }

                    break;
                }
            }

            // Execute propagation step
            //
            // Loop over all directions
#pragma unroll
            for (int dir = 0; dir < this->NUM_DIR; dir++)
            {
                // Reset the memory offset.
                int offset = 0;

                // The cell is located in a row with even index value
                if (pos_y % 2 == 0)
                {
                    // Construct the correct memory offset
                    //
                    // Apply a default offset value
                    offset += m_model->offset_to_neighbor_even[dir];

                    // Correct the offset in the current direction if the cell is located on boundaries
                    if (on_eastern_boundary) {

                        offset += m_model->offset_to_western_boundary_even[dir];
                    }

                    if (on_northern_boundary) {

                        offset += m_model->offset_to_southern_boundary_even[dir];
                    }

                    if (on_western_boundary) {

                        offset += m_model->offset_to_eastern_boundary_even[dir];
                    }

                    if (on_southern_boundary) {

                        offset += m_model->offset_to_northern_boundary_even[dir];
                    }

                // The cell is located in a row with odd index value
                } else if (pos_y % 2 != 0) {

                    // Construct the correct memory offset.
                    //
                    // Apply a default offset value
                    offset += m_model->offset_to_neighbor_odd[dir];

                    // Correct the offset in the current direction if the cell is located on boundaries
                    if (on_eastern_boundary) {

                        offset += m_model->offset_to_western_boundary_odd[dir];
                    }

                    if (on_northern_boundary) {

                        offset += m_model->offset_to_southern_boundary_odd[dir];
                    }

                    if (on_western_boundary) {

                        offset += m_model->offset_to_eastern_boundary_odd[dir];
                    }

                    if (on_southern_boundary) {

                        offset += m_model->offset_to_northern_boundary_odd[dir];
                    }
                }

                // Push the states of the cell to its "neighbor" cells in the different directions
                m_node_state_tmp_cpu[node_idx[dir] + offset * this->NUM_DIR] = bool(node_state_tmp[dir]);
            }

        } /* FOR cell */

    }}); /* FOR block */

    // Update the node states.
    auto node_state_cpu_tmp = this->m_node_state_cpu.ptr();
    this->m_node_state_cpu = m_node_state_tmp_cpu.ptr();
    m_node_state_tmp_cpu = node_state_cpu_tmp;
}

// Applies a body force in the specified direction (x or y) and with the
// specified intensity to the particles. E.g., if the intensity is equal 100,
// every 100th particle changes it's direction, if feasible.
template<Model model_>
void OMP_Lattice<model_>:: OMP_Lattice::apply_body_force(const int forcing) {

    // TODO: Set the seed for the random number generation on the device.
    int seed = time(NULL);

    // Set a maximum number of iterations to find particles which can be reverted.
    const int it_max = 2 * this->m_num_cells;

    // Set the number of iterations to zero.
    int it = 0;

    // Number of particles which have been reverted.
    int reverted_particles = 0;

    // Loop over all cells.
    do
    {
        int cell = rand() % this->m_num_cells;

    	it++;

        // Get the type of the cell, i.e. fluid or solid.
        // Note that body forces are applied to fluid cells only.
        char cell_type = this->m_cell_type_cpu[cell];

        // Check weather the cell working on is a fluid cell.
        if (cell_type == 0) {

            // Define an array for the global indices of the nodes in the cell.
            int node_idx[this->NUM_DIR];

            // Define an array for the states of the nodes in the cell.
            char node_state[this->NUM_DIR];

            // The thread working on the cell has to know about the states of the
            // nodes within the cell, therefore looping over all directions and
            // look it up.
#pragma unroll
            for (int dir = 0; dir < this->NUM_DIR; ++dir) {

                node_idx[dir] = dir + cell * this->NUM_DIR;
                node_state[dir] = bool(this->m_node_state_cpu[node_idx[dir]]);
            }

            // Create a temporary array to copy the node states into.
            char node_state_tmp[this->NUM_DIR];

            // Copy the current states of the nodes to the temporary array.
#pragma unroll
            for (int dir = 0; dir < this->NUM_DIR; ++dir) {

                node_state_tmp[dir] = node_state[dir];
            }

            if (model_ == Model::HPP) {

                if (this->m_bf_dir == 'x' && (node_state[0] == 0) && (node_state[2] == 1)) {

					node_state_tmp[0] = 1;
					node_state_tmp[2] = 0;

					reverted_particles++;

                } else if (this->m_bf_dir == 'y' && (node_state[1] == 1) && (node_state[3] == 0)) {

					node_state_tmp[1] = 0;
					node_state_tmp[3] = 1;

					reverted_particles++;
				}
        	}

            else if (model_ == Model::FHP) {

                if (this->m_bf_dir == 'x' && (node_state[0] == 0) && (node_state[3] == 1)) {

					node_state_tmp[0] = 1;
					node_state_tmp[3] = 0;

					reverted_particles++;

                } else if (this->m_bf_dir == 'y') {

					if ((node_state[1] == 1) && (node_state[5] == 0)) {

						node_state_tmp[1] = 0;
						node_state_tmp[5] = 1;

						reverted_particles++;
					}

					if ((node_state[2] == 1) && (node_state[4] == 0)) {

						node_state_tmp[2] = 0;
						node_state_tmp[4] = 1;

						reverted_particles++;
					}
				}
        	}

//            // Loop over all directions.
//#pragma unroll
//            for (int dir = 0; dir < n_dir; ++dir) {
//
//                // Body force acting in x direction.
//                if (bf_dir == 'x') {
//
//					// TODO: Exchange the states of the nodes with the the states of
//					//       the mirrored directions along the y axis if feasible.
//					if ((fabs(lattice_vec_x[dir] - 1.0) < 1.0e-06) &&
//						(node_state[dir] < node_state[mirrored_dir_y[dir]])) {
//
//						node_state_tmp[dir                ] = node_state[mirrored_dir_y[dir]];
//						node_state_tmp[mirrored_dir_y[dir]] = node_state[dir                ];
//					}
//                }
//
//                // Body force acting in y direction.
//                else if (bf_dir == 'y') {
//
//					// TODO: Exchange the states of the nodes with the the states of
//					//       the mirrored directions along the x axis if feasible.
//					if ((lattice_vec_y[dir] < 1.0e-06) &&
//						(node_state[dir] < node_state[mirrored_dir_x[dir]])) {
//
//						node_state_tmp[dir                ] = node_state[mirrored_dir_x[dir]];
//						node_state_tmp[mirrored_dir_x[dir]] = node_state[dir                ];
//					}
//                }
//
//#ifdef DEBUG
//                // Invalid body force direction.
//                else {
//
//                    printf("ERROR in apply_body_force(): "
//                           "Invalid body force direction %c.\n", bf_dir);
//                }
//#endif
//            }

            // Write the new node states back to the data array.
            //
            // Loop over all directions.
#pragma unroll
            for (int dir = 0; dir < this->NUM_DIR; dir++)
            {
                this->m_node_state_cpu[node_idx[dir]] = bool(node_state_tmp[dir]);
            }

        } /* IF cell_type */

    } while ((reverted_particles < forcing) && (it < it_max));
}

// Computes quantities of interest as a post-processing procedure.
template<Model model_>
void OMP_Lattice<model_>::post_process() {

	// Computes cell quantities of interest as a post-processing procedure.
	cell_post_process();

	// Computes coarse grained quantities of interest as a post-processing procedure.
	mean_post_process();
}

// Computes cell quantities of interest as a post-processing procedure.
template<Model model_>
void OMP_Lattice<model_>::cell_post_process()
{
    // Loop over lattice cells
    tbb::parallel_for(tbb::blocked_range<int>(0, this->m_num_cells), [&](const tbb::blocked_range<int>& r) {
    for (int cell = r.begin(); cell != r.end(); ++cell)
    {
        // Initialize the cell quantities to be computed
        char cell_density    = 0;
        Real cell_momentum_x = 0.0;
        Real cell_momentum_y = 0.0;

        // Loop over nodes within the current cell
#pragma unroll
        for (int dir = 0; dir < this->NUM_DIR; ++dir) {

            char node_state = bool(this->m_node_state_out_cpu[dir + cell * this->NUM_DIR]);

            // Sum up the node states
            cell_density += node_state;

            // Sum up the node states multiplied by the lattice vector component for the current
            // direction
            cell_momentum_x += node_state * ModelDesc::LATTICE_VEC_X[dir];
            cell_momentum_y += node_state * ModelDesc::LATTICE_VEC_Y[dir];
        }

        // Write the computed cell quantities to the related data arrays
        this->m_cell_density_cpu [cell                        ] = (Real) cell_density;
        this->m_cell_momentum_cpu[cell * this->SPATIAL_DIM    ] = cell_momentum_x;
        this->m_cell_momentum_cpu[cell * this->SPATIAL_DIM + 1] = cell_momentum_y;

    } // for cell
    });
}

// Computes coarse grained quantities of interest as a post-processing procedure
template<Model model_>
void OMP_Lattice<model_>::mean_post_process()
{
    const int r = this->m_coarse_graining_radius;

    tbb::parallel_for(tbb::blocked_range<int>(0, this->m_num_coarse_cells), [&](const tbb::blocked_range<int>& range) {
    for (int coarse_cell = range.begin(); coarse_cell != range.end(); ++coarse_cell)
    {
        // Get cell in the middle of the coarse cell
        const int cell = r + r * this->m_dim_x
                + (coarse_cell % this->m_coarse_dim_x) * (2 * r + 1)
                + (coarse_cell / this->m_coarse_dim_x) * (2 * r + 1) * this->m_dim_x;

        // Calculate the position of the cell in x direction
        int pos_x = cell % this->m_dim_x;

        // Initialize the coarse grained quantities to be computed
        Real mean_density    = 0.0;
        Real mean_momentum_x = 0.0;
        Real mean_momentum_y = 0.0;

        // Initialize the number of actual existing coarse graining neighbor cells
        int n_exist_neighbors = 0;

        // The thread working on the cell has to know the cell quantities of the coarse graining
        // neighbor cells, therefore looping over all neighbor cells and look it up
#pragma unroll
        for (int y = -r; y <= r; ++y) {

            for (int x = -r; x <= r; ++x) {

                // Get the index of the coarse graining neighbor cell
                int neighbor_idx = cell + y * this->m_dim_x + x;

                // Get the position of the coarse graining neighbor cell in x direction
                int pos_x_neighbor = neighbor_idx % this->m_dim_x;

                // Check weather the coarse graining neighbor cell is valid
                if ((neighbor_idx >= 0) &&
                    (neighbor_idx < this->m_num_cells) &&
                    (abs(pos_x_neighbor - pos_x) <= r)) {

                    // Increase the number of existing coarse graining neighbor cells
                    n_exist_neighbors++;

                    mean_density    += this->m_cell_density_cpu [neighbor_idx                        ];
                    mean_momentum_x += this->m_cell_momentum_cpu[neighbor_idx * this->SPATIAL_DIM    ];
                    mean_momentum_y += this->m_cell_momentum_cpu[neighbor_idx * this->SPATIAL_DIM + 1];
                }
            }
        }

        // Write the computed coarse grained quantities to the related data arrays
        this->m_mean_density_cpu [coarse_cell                        ] = mean_density    / ((Real) n_exist_neighbors);
        this->m_mean_momentum_cpu[coarse_cell * this->SPATIAL_DIM    ] = mean_momentum_x / ((Real) n_exist_neighbors);
        this->m_mean_momentum_cpu[coarse_cell * this->SPATIAL_DIM + 1] = mean_momentum_y / ((Real) n_exist_neighbors);

    }}); // for coarse_cell
}

// Allocates the memory for the arrays on the host (CPU)
template<Model model_>
void OMP_Lattice<model_>::allocate_memory()
{
    // Allocate host memory
    cu_verify(cudaMallocHost((void **) &this->m_cell_type_cpu,                           this->m_num_cells        * sizeof(char)));
    cu_verify(cudaMallocHost((void **) &this->m_cell_density_cpu,                        this->m_num_cells        * sizeof(Real)));
    cu_verify(cudaMallocHost((void **) &this->m_mean_density_cpu,                        this->m_num_coarse_cells * sizeof(Real)));
    cu_verify(cudaMallocHost((void **) &this->m_cell_momentum_cpu,   this->SPATIAL_DIM * this->m_num_cells        * sizeof(Real)));
    cu_verify(cudaMallocHost((void **) &this->m_mean_momentum_cpu,   this->SPATIAL_DIM * this->m_num_coarse_cells * sizeof(Real)));

    this->m_node_state_cpu.resize    (this->m_num_nodes);
          m_node_state_tmp_cpu.resize(this->m_num_nodes);
    this->m_node_state_out_cpu.resize(this->m_num_nodes);
}

// Frees the memory for the arrays on the host (CPU)
template<Model model_>
void OMP_Lattice<model_>::free_memory()
{
    // Free CPU memory
    cu_verify(cudaFreeHost(this->m_cell_type_cpu));
    cu_verify(cudaFreeHost(this->m_cell_density_cpu));
    cu_verify(cudaFreeHost(this->m_mean_density_cpu));
    cu_verify(cudaFreeHost(this->m_cell_momentum_cpu));
    cu_verify(cudaFreeHost(this->m_mean_momentum_cpu));

    this->m_cell_type_cpu       = NULL;
    this->m_cell_density_cpu    = NULL;
    this->m_mean_density_cpu    = NULL;
    this->m_cell_momentum_cpu   = NULL;
    this->m_mean_momentum_cpu   = NULL;
}

// Sets (proper) parallelization parameters
template<Model model_>
void OMP_Lattice<model_>::setup_parallel()
{

#pragma omp parallel for
	for (unsigned int i = 0; i < omp_get_num_threads(); ++i)
	{
		if (omp_get_thread_num() == 0)
		{
			// Only execute in master thread.
            printf("OMP configuration parameters: Executing calculation with %d threads.\n\n", omp_get_num_threads());
		}
	}
}

// Computes the mean velocity of the lattice.
template<Model model_>
std::vector<Real> OMP_Lattice<model_>::get_mean_velocity() {

    std::vector<Real> mean_velocity(this->SPATIAL_DIM, 0.0);

    Real sum_x_vel = 0.0;
    Real sum_y_vel = 0.0;

    unsigned int counter = 0;

    // Sum up all (fluid) cell x and y velocity components.
#pragma omp parallel for reduction(+: sum_x_vel, sum_y_vel)
    for (unsigned int n = 0; n < this->m_num_cells; ++n) {

        if (this->m_cell_type_cpu[n] == 0) {

        	counter++;

            Real cell_density = this->m_cell_density_cpu[n];

			if (cell_density > 1.0e-06) {

                sum_x_vel += this->m_cell_momentum_cpu[n * this->SPATIAL_DIM    ] / cell_density;
                sum_y_vel += this->m_cell_momentum_cpu[n * this->SPATIAL_DIM + 1] / cell_density;
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
template class OMP_Lattice<Model::FHP>;

} // namespace lgca
