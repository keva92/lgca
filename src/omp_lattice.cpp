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

// Creates a CUDA parallelized lattice gas cellular automaton object
// of the specified properties.
OMP_Lattice::OMP_Lattice(const string test_case,
                         const Real Re, const Real Ma_s,
                         const int n_dir,
                         const int coarse_graining_radius)

			             : Lattice(test_case, Re, Ma_s, n_dir, coarse_graining_radius) {

    // Allocate the memory for the arrays on the host (CPU).
    allocate_memory();

	// Memory offset to neighbor cells in the different directions for the
	// propagation step.
	offset_to_neighbor_even = new int[n_dir];
	offset_to_neighbor_odd =  new int[n_dir];

	// Memory offset to related cells of the opposite boundary in the different
	// directions in case of periodic boundaries.
	offset_to_eastern_boundary_even  = new int[n_dir];
	offset_to_eastern_boundary_odd   = new int[n_dir];
	offset_to_northern_boundary_even = new int[n_dir];
	offset_to_northern_boundary_odd  = new int[n_dir];
	offset_to_western_boundary_even  = new int[n_dir];
	offset_to_western_boundary_odd   = new int[n_dir];
	offset_to_southern_boundary_even = new int[n_dir];
	offset_to_southern_boundary_odd  = new int[n_dir];

	// Inverse direction indices for each lattice direction.
	inverse_dir = new char[n_dir];

	// Mirrored direction indices for each lattice direction with respect
	// to the x and y axis.
	mirrored_dir_x = new char[n_dir];
	mirrored_dir_y = new char[n_dir];

	// Set the components of the lattice vectors for the different directions.
	//
	// Set the model based values according to the number of lattice directions.
	switch (n_dir) {

		// HPP model.
		case 4:
		{
			// Define the memory offsets in the different directions for
			// cells in rows with even and odd indices (which are the same
			// in case of HPP).
			//
			// The cell is located in a row with even index value.
			offset_to_neighbor_even[0] = 1;
			offset_to_neighbor_even[1] = n_x;
			offset_to_neighbor_even[2] = -1;
			offset_to_neighbor_even[3] = -n_x;

			offset_to_eastern_boundary_even[0] = 0;
			offset_to_eastern_boundary_even[1] = 0;
			offset_to_eastern_boundary_even[2] = n_x;
			offset_to_eastern_boundary_even[3] = 0;

			offset_to_northern_boundary_even[0] = 0;
			offset_to_northern_boundary_even[1] = 0;
			offset_to_northern_boundary_even[2] = 0;
			offset_to_northern_boundary_even[3] = n_x * n_y;

			offset_to_western_boundary_even[0] = -n_x;
			offset_to_western_boundary_even[1] = 0;
			offset_to_western_boundary_even[2] = 0;
			offset_to_western_boundary_even[3] = 0;

			offset_to_southern_boundary_even[0] = 0;
			offset_to_southern_boundary_even[1] = -n_x * n_y;
			offset_to_southern_boundary_even[2] = 0;
			offset_to_southern_boundary_even[3] = 0;

			// The cell is located in a row with odd index value.
			offset_to_neighbor_odd[0] = offset_to_neighbor_even[0];
			offset_to_neighbor_odd[1] = offset_to_neighbor_even[1];
			offset_to_neighbor_odd[2] = offset_to_neighbor_even[2];
			offset_to_neighbor_odd[3] = offset_to_neighbor_even[3];

			offset_to_eastern_boundary_odd[0] = offset_to_eastern_boundary_even[0];
			offset_to_eastern_boundary_odd[1] = offset_to_eastern_boundary_even[1];
			offset_to_eastern_boundary_odd[2] = offset_to_eastern_boundary_even[2];
			offset_to_eastern_boundary_odd[3] = offset_to_eastern_boundary_even[3];

			offset_to_northern_boundary_odd[0] = offset_to_northern_boundary_even[0];
			offset_to_northern_boundary_odd[1] = offset_to_northern_boundary_even[1];
			offset_to_northern_boundary_odd[2] = offset_to_northern_boundary_even[2];
			offset_to_northern_boundary_odd[3] = offset_to_northern_boundary_even[3];

			offset_to_western_boundary_odd[0] = offset_to_western_boundary_even[0];
			offset_to_western_boundary_odd[1] = offset_to_western_boundary_even[1];
			offset_to_western_boundary_odd[2] = offset_to_western_boundary_even[2];
			offset_to_western_boundary_odd[3] = offset_to_western_boundary_even[3];

			offset_to_southern_boundary_odd[0] = offset_to_southern_boundary_even[0];
			offset_to_southern_boundary_odd[1] = offset_to_southern_boundary_even[1];
			offset_to_southern_boundary_odd[2] = offset_to_southern_boundary_even[2];
			offset_to_southern_boundary_odd[3] = offset_to_southern_boundary_even[3];

			inverse_dir[0] = 2;
			inverse_dir[1] = 3;
			inverse_dir[2] = 0;
			inverse_dir[3] = 1;

			mirrored_dir_x[0] = 0;
			mirrored_dir_x[1] = 3;
			mirrored_dir_x[2] = 2;
			mirrored_dir_x[3] = 1;

			mirrored_dir_y[0] = 2;
			mirrored_dir_y[1] = 1;
			mirrored_dir_y[2] = 0;
			mirrored_dir_y[3] = 3;

			break;
		}

		// FHP model.
		case 6:
		{
			// Define the memory offsets in the different directions for
			// cells in rows with even and odd indices.
			//
			// The cell is located in a row with even index value.
			offset_to_neighbor_even[0] = 1;
			offset_to_neighbor_even[1] = n_x;
			offset_to_neighbor_even[2] = n_x - 1;
			offset_to_neighbor_even[3] = -1;
			offset_to_neighbor_even[4] = -n_x - 1;
			offset_to_neighbor_even[5] = -n_x;

			offset_to_eastern_boundary_even[0] = 0;
			offset_to_eastern_boundary_even[1] = 0;
			offset_to_eastern_boundary_even[2] = n_x;
			offset_to_eastern_boundary_even[3] = n_x;
			offset_to_eastern_boundary_even[4] = n_x;
			offset_to_eastern_boundary_even[5] = 0;

			offset_to_northern_boundary_even[0] = 0;
			offset_to_northern_boundary_even[1] = 0;
			offset_to_northern_boundary_even[2] = 0;
			offset_to_northern_boundary_even[3] = 0;
			offset_to_northern_boundary_even[4] = n_x * n_y;
			offset_to_northern_boundary_even[5] = n_x * n_y;

			offset_to_western_boundary_even[0] = -n_x;
			offset_to_western_boundary_even[1] = 0;
			offset_to_western_boundary_even[2] = 0;
			offset_to_western_boundary_even[3] = 0;
			offset_to_western_boundary_even[4] = 0;
			offset_to_western_boundary_even[5] = 0;

			offset_to_southern_boundary_even[0] = 0;
			offset_to_southern_boundary_even[1] = -n_x * n_y;
			offset_to_southern_boundary_even[2] = -n_x * n_y + 1;
			offset_to_southern_boundary_even[3] = 0;
			offset_to_southern_boundary_even[4] = 0;
			offset_to_southern_boundary_even[5] = 0;

			// The cell is located in a row with odd index value.
			offset_to_neighbor_odd[0] = 1;
			offset_to_neighbor_odd[1] = n_x + 1;
			offset_to_neighbor_odd[2] = n_x;
			offset_to_neighbor_odd[3] = -1;
			offset_to_neighbor_odd[4] = -n_x;
			offset_to_neighbor_odd[5] = -n_x + 1;

			offset_to_eastern_boundary_odd[0] = 0;
			offset_to_eastern_boundary_odd[1] = 0;
			offset_to_eastern_boundary_odd[2] = 0;
			offset_to_eastern_boundary_odd[3] = n_x;
			offset_to_eastern_boundary_odd[4] = 0;
			offset_to_eastern_boundary_odd[5] = 0;

			offset_to_northern_boundary_odd[0] = 0;
			offset_to_northern_boundary_odd[1] = 0;
			offset_to_northern_boundary_odd[2] = 0;
			offset_to_northern_boundary_odd[3] = 0;
			offset_to_northern_boundary_odd[4] = n_x * n_y;
			offset_to_northern_boundary_odd[5] = n_x * n_y;

			offset_to_western_boundary_odd[0] = -n_x;
			offset_to_western_boundary_odd[1] = -n_x;
			offset_to_western_boundary_odd[2] = 0;
			offset_to_western_boundary_odd[3] = 0;
			offset_to_western_boundary_odd[4] = 0;
			offset_to_western_boundary_odd[5] = -n_x;

			offset_to_southern_boundary_odd[0] = 0;
			offset_to_southern_boundary_odd[1] = -n_x * n_y;
			offset_to_southern_boundary_odd[2] = -n_x * n_y;
			offset_to_southern_boundary_odd[3] = 0;
			offset_to_southern_boundary_odd[4] = 0;
			offset_to_southern_boundary_odd[5] = 0;

			inverse_dir[0] = 3;
			inverse_dir[1] = 4;
			inverse_dir[2] = 5;
			inverse_dir[3] = 0;
			inverse_dir[4] = 1;
			inverse_dir[5] = 2;

			mirrored_dir_x[0] = 0;
			mirrored_dir_x[1] = 5;
			mirrored_dir_x[2] = 4;
			mirrored_dir_x[3] = 3;
			mirrored_dir_x[4] = 2;
			mirrored_dir_x[5] = 1;

			mirrored_dir_y[0] = 3;
			mirrored_dir_y[1] = 2;
			mirrored_dir_y[2] = 1;
			mirrored_dir_y[3] = 0;
			mirrored_dir_y[4] = 5;
			mirrored_dir_y[5] = 4;

			break;
		}

#ifdef DEBUG
		default:
		{
			printf("ERROR in OMP_Lattice::OMP_Lattice: Invalid number of directions %d!\n", n_dir);
			abort();
			break;
		}
#endif

	}

	// Lattice vector components in the different directions.
    lattice_vec_x = new Real[n_dir];
    lattice_vec_y = new Real[n_dir];

    // Set the components of the lattice vectors for the different directions.
    //
    // Loop over all directions.
    for (int dir = 0; dir < n_dir; ++dir) {

        lattice_vec_x[dir] = cos(2.0 * M_PI / ((Real) n_dir) * ((Real) dir));
        lattice_vec_y[dir] = sin(2.0 * M_PI / ((Real) n_dir) * ((Real) dir));
    }
}

// Deletes the openMP parallelized lattice gas cellular automaton object.
OMP_Lattice::~OMP_Lattice() {

	free_memory();
}

// Performs the collision and propagation step on the lattice gas automaton.
void OMP_Lattice::collide_and_propagate(unsigned int step) {

    // TODO: Set the seed for the random number generation.
    int seed = time(NULL);

#ifdef DEBUG
            // Check weather the domain dimensions are valid for the FHP model.
            if (n_y % 2 != 0 && n_dir == 6) {

                printf("ERROR in collide_and_propagate_kernel(): "
                       "Invalid domain dimension in y direction.\n");
                abort();
            }
#endif

    // Loop over all cells.
//#pragma omp parallel for
//	for (unsigned int cell = 0; cell < n_cells; ++cell)
//	{
    tbb::parallel_for(tbb::blocked_range<int>(0, n_cells), [&](const tbb::blocked_range<int>& r) {
    for (int cell = r.begin(); cell != r.end(); ++cell)
    {
		// Calculate the position of the cell in x direction (column index).
		int pos_x = cell % n_x;

		// Calculate the position of the cell in y direction (row index).
		int pos_y = cell / n_x;

		// Get the type of the cell, i.e. fluid or solid.
		// This has to be taken into account during the collision step, where
		// cells behave different according to their type.
		char cell_type = cell_type_cpu[cell];

		// Check weather the cell is located on boundaries.
		bool on_eastern_boundary  = (cell + 1) % n_x == 0;
		bool on_northern_boundary = cell >= (n_cells - n_x);
		bool on_western_boundary  = cell % n_x == 0;
		bool on_southern_boundary = cell < n_x;

		// Define an array for the global indices of the nodes in the cell.
		int node_idx[n_dir];

		// Define an array for the states of the nodes in the cell.
		char node_state[n_dir];

		// Execute collision step.
		//
		// The thread working on the cell has to know about the states of the
		// nodes within the cell, therefore looping over all directions and
		// look it up.
	#pragma unroll
		for (int dir = 0; dir < n_dir; ++dir) {

			node_idx[dir] = cell + dir * n_cells;
			node_state[dir] = node_state_cpu[node_idx[dir]];
		}

		// TODO: Create a random boolean value for the collision step.
		// bool rand_bool = cu_random_bool(seed, cell);
		bool rand_bool =       ((pos_x % 2) == (pos_y % 2))
						 - 1 * ((pos_x % 2) == (pos_y % 2)) * (step % 2)
						 + 1 * ((pos_x % 2) != (pos_y % 2)) * (step % 2);

		// Create a temporary array to copy the node states into.
		char node_state_tmp[n_dir];

		// Copy the actual states of the nodes to the temporary array.
	#pragma unroll
		for (int dir = 0; dir < n_dir; ++dir) {

			node_state_tmp[dir] = node_state[dir];
		}

		switch (cell_type) {

			// The cell working on is a fluid cell ("normal" collision).
			case 0:
			{
				// Using the the HPP model.
				if (n_dir == 4) {

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
				} else if (n_dir == 6) {

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

	#ifdef DEBUG
				else {

					printf("ERROR in OMP_Lattice::collide_and_propagate(): "
						   "Invalid number of directions %d.\n", n_dir);
				}
	#endif

				break;
			}

			// The cell working on is a solid cell of bounce back type.
			case 1:
			{
				// Loop over all directions.
				// #pragma unroll
				for (int dir = 0; dir < n_dir; ++dir) {

					// Exchange the states of the nodes with the the states of
					// the inverse directions.
					node_state_tmp[dir] = node_state[inverse_dir[dir]];
				}

				break;
			}

			// TODO: The cell working on is a solid cell of bounce forward type.
			case 2:
			{
				// Loop over all directions.
	#pragma unroll
				for (int dir = 0; dir < n_dir; ++dir) {

					if (on_northern_boundary || on_southern_boundary) {

						// Exchange the states of the nodes with the the states of
						// the mirrored directions along the x axis.
						node_state_tmp[dir] = node_state[mirrored_dir_x[dir]];
					}

					if (on_eastern_boundary || on_western_boundary) {

						// Exchange the states of the nodes with the the states of
						// the mirrored directions along the y axis.
						node_state_tmp[dir] = node_state[mirrored_dir_y[dir]];
					}
				}

				break;
			}

	#ifdef DEBUG
			// Invalid cell type.
			default:
			{
				printf("ERROR in OMP_Lattice::collide_and_propagate(): "
					   "Invalid cell type %d.\n", cell_type);
				break;
			}
	#endif

		}

		// Execute propagation step.
		//
		// Loop over all directions.
#pragma unroll
		for (int dir = 0; dir < n_dir; dir++)
		{
			// Reset the memory offset.
			int offset = 0;

            // The cell is located in a row with even index value.
            if (pos_y % 2 == 0)
            {
				// Construct the correct memory offset.
				//
				// Apply a default offset value.
				offset += offset_to_neighbor_even[dir];

				// Correct the offset in the current direction if the cell is
				// located on boundaries.
				if (on_eastern_boundary) {

					offset += offset_to_western_boundary_even[dir];
				}

				if (on_northern_boundary) {

					offset += offset_to_southern_boundary_even[dir];
				}

				if (on_western_boundary) {

					offset += offset_to_eastern_boundary_even[dir];
				}

				if (on_southern_boundary) {

					offset += offset_to_northern_boundary_even[dir];
				}

            // The cell is located in a row with odd index value.
            } else if (pos_y % 2 != 0) {

				// Construct the correct memory offset.
				//
				// Apply a default offset value.
				offset += offset_to_neighbor_odd[dir];

				// Correct the offset in the current direction if the cell is
				// located on boundaries.
				if (on_eastern_boundary) {

					offset += offset_to_western_boundary_odd[dir];
				}

				if (on_northern_boundary) {

					offset += offset_to_southern_boundary_odd[dir];
				}

				if (on_western_boundary) {

					offset += offset_to_eastern_boundary_odd[dir];
				}

				if (on_southern_boundary) {

					offset += offset_to_northern_boundary_odd[dir];
				}
            }

			// Push the states of the cell to its "neighbor" cells in the
			// different directions.
			node_state_tmp_cpu[node_idx[dir] + offset] = node_state_tmp[dir];
		}

    } /* FOR cell */
    });

    // Update the node states.
    char* node_state_cpu_tmp = node_state_cpu;
    node_state_cpu = node_state_tmp_cpu;
    node_state_tmp_cpu = node_state_cpu_tmp;
}

// Applies a body force in the specified direction (x or y) and with the
// specified intensity to the particles. E.g., if the intensity is equal 100,
// every 100th particle changes it's direction, if feasible.
void OMP_Lattice::apply_body_force(const int forcing) {

    // TODO: Set the seed for the random number generation on the device.
    int seed = time(NULL);

    // Set a maximum number of iterations to find particles which can be reverted.
    const int it_max = 2 * n_cells;

    // Set the number of iterations to zero.
    int it = 0;

    // Number of particles which have been reverted.
    int reverted_particles = 0;

    // Loop over all cells.
    do
    {
    	int cell = rand() % n_cells;

    	it++;

        // Get the type of the cell, i.e. fluid or solid.
        // Note that body forces are applied to fluid cells only.
        char cell_type = cell_type_cpu[cell];

        // Check weather the cell working on is a fluid cell.
        if (cell_type == 0) {

            // Define an array for the global indices of the nodes in the cell.
            int node_idx[n_dir];

            // Define an array for the states of the nodes in the cell.
            char node_state[n_dir];

            // The thread working on the cell has to know about the states of the
            // nodes within the cell, therefore looping over all directions and
            // look it up.
    #pragma unroll
            for (int dir = 0; dir < n_dir; ++dir) {

                node_idx[dir] = cell + dir * n_cells;
                node_state[dir] = node_state_cpu[node_idx[dir]];
            }

            // Create a temporary array to copy the node states into.
            char node_state_tmp[n_dir];

            // Copy the current states of the nodes to the temporary array.
    #pragma unroll
            for (int dir = 0; dir < n_dir; ++dir) {

                node_state_tmp[dir] = node_state[dir];
            }

        	if (n_dir == 4) {

				if (bf_dir == 'x' && (node_state[0] == 0) && (node_state[2] == 1)) {

					node_state_tmp[0] = 1;
					node_state_tmp[2] = 0;

					reverted_particles++;

				} else if (bf_dir == 'y' && (node_state[1] == 1) && (node_state[3] == 0)) {

					node_state_tmp[1] = 0;
					node_state_tmp[3] = 1;

					reverted_particles++;
				}
        	}

        	else if (n_dir == 6) {

				if (bf_dir == 'x' && (node_state[0] == 0) && (node_state[3] == 1)) {

					node_state_tmp[0] = 1;
					node_state_tmp[3] = 0;

					reverted_particles++;

				} else if (bf_dir == 'y') {

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
            for(int dir = 0; dir < n_dir; dir++)
            {
                node_state_cpu[node_idx[dir]] = node_state_tmp[dir];
            }

        } /* IF cell_type */

    } while ((reverted_particles < forcing) && (it < it_max));
}

// Computes quantities of interest as a post-processing procedure.
void OMP_Lattice::post_process() {

	// Computes cell quantities of interest as a post-processing procedure.
	cell_post_process();

	// Computes coarse grained quantities of interest as a post-processing procedure.
	mean_post_process();
}

// Computes cell quantities of interest as a post-processing procedure.
void OMP_Lattice::cell_post_process()
{
    // Loop over lattice cells
//#pragma omp parallel for
//    for (int cell = 0; cell < n_cells; ++cell)
//    {
    tbb::parallel_for(tbb::blocked_range<int>(0, n_cells), [&](const tbb::blocked_range<int>& r) {
    for (int cell = r.begin(); cell != r.end(); ++cell)
    {
        // Initialize the cell quantities to be computed
        char cell_density    = 0;
        Real cell_momentum_x = 0.0;
        Real cell_momentum_y = 0.0;

        // Loop over nodes within the current cell
#pragma unroll
        for (int dir = 0; dir < n_dir; ++dir) {

            int node_idx = cell + dir * n_cells;
            char node_state = node_state_out_cpu[node_idx];

            // Sum up the node states
            cell_density += node_state;

            // Sum up the node states multiplied by the lattice vector component for the current
            // direction
            cell_momentum_x += node_state * lattice_vec_x[dir];
            cell_momentum_y += node_state * lattice_vec_y[dir];
        }

        // Write the computed cell quantities to the related data arrays
        cell_density_cpu [cell          ] = (Real) cell_density;
        cell_momentum_cpu[cell          ] =        cell_momentum_x;
        cell_momentum_cpu[cell + n_cells] =        cell_momentum_y;

    } // for cell
    });
}

// Computes coarse grained quantities of interest as a post-processing procedure
void OMP_Lattice::mean_post_process()
{
    // Loop over all cells
//#pragma omp parallel for
//    for (unsigned int cell = 0; cell < n_cells; ++cell)
//    {
    tbb::parallel_for(tbb::blocked_range<int>(0, n_cells), [&](const tbb::blocked_range<int>& r) {
    for (int cell = r.begin(); cell != r.end(); ++cell)
    {
        // Calculate the position of the cell in x direction
        int pos_x = cell % n_x;

        // Initialize the coarse grained quantities to be computed
        Real mean_density    = 0.0;
        Real mean_momentum_x = 0.0;
        Real mean_momentum_y = 0.0;

        // Initialize the number of actual existing coarse graining neighbor cells
        int n_exist_neighbors = 0;

        // The thread working on the cell has to know the cell quantities of the coarse graining
        // neighbor cells, therefore looping over all neighbor cells and look it up
#pragma unroll
        for (int y = -coarse_graining_radius; y <= coarse_graining_radius; ++y) {

            for (int x = -coarse_graining_radius; x <= coarse_graining_radius; ++x) {

                // Get the index of the coarse graining neighbor cell
                int neighbor_idx = cell + y * n_x + x;

                // Get the position of the coarse graining neighbor cell in x direction
                int pos_x_neighbor = neighbor_idx % n_x;

                // Check weather the coarse graining neighbor cell is valid
                if ((neighbor_idx >= 0) &
                    (neighbor_idx < n_cells) &&
                    (abs(pos_x_neighbor - pos_x) <= coarse_graining_radius)) {

                    // Increase the number of existing coarse graining neighbor cells
                    n_exist_neighbors++;

                    mean_density    += cell_density_cpu [neighbor_idx          ];
                    mean_momentum_x += cell_momentum_cpu[neighbor_idx          ];
                    mean_momentum_y += cell_momentum_cpu[neighbor_idx + n_cells];
                }
            }
        }

        // Write the computed coarse grained quantities to the related data arrays
        mean_density_cpu [cell          ] = mean_density    / ((Real) n_exist_neighbors);
        mean_momentum_cpu[cell          ] = mean_momentum_x / ((Real) n_exist_neighbors);
        mean_momentum_cpu[cell + n_cells] = mean_momentum_y / ((Real) n_exist_neighbors);

    } // for cell
    });
}

// Allocates the memory for the arrays on the host (CPU).
void OMP_Lattice::allocate_memory()
{
    // Allocate host memory.
    cu_verify(cudaMallocHost((void **) &node_state_cpu,          n_nodes * sizeof(char)));
    cu_verify(cudaMallocHost((void **) &node_state_tmp_cpu,      n_nodes * sizeof(char)));
    cu_verify(cudaMallocHost((void **) &node_state_out_cpu,      n_nodes * sizeof(char)));
    cu_verify(cudaMallocHost((void **) &cell_type_cpu,           n_cells * sizeof(char)));
    cu_verify(cudaMallocHost((void **) &cell_density_cpu,        n_cells * sizeof(Real)));
    cu_verify(cudaMallocHost((void **) &mean_density_cpu,        n_cells * sizeof(Real)));
    cu_verify(cudaMallocHost((void **) &cell_momentum_cpu, dim * n_cells * sizeof(Real)));
    cu_verify(cudaMallocHost((void **) &mean_momentum_cpu, dim * n_cells * sizeof(Real)));
}

// Frees the memory for the arrays on the host (CPU).
void OMP_Lattice::free_memory()
{
    // Free CPU memory.
    cu_verify(cudaFreeHost(node_state_cpu));
    cu_verify(cudaFreeHost(node_state_tmp_cpu));
    cu_verify(cudaFreeHost(node_state_out_cpu));
    cu_verify(cudaFreeHost(cell_type_cpu));
    cu_verify(cudaFreeHost(cell_density_cpu));
    cu_verify(cudaFreeHost(mean_density_cpu));
    cu_verify(cudaFreeHost(cell_momentum_cpu));
    cu_verify(cudaFreeHost(mean_momentum_cpu));

	delete[] offset_to_neighbor_even;
	delete[] offset_to_neighbor_odd;
	delete[] offset_to_eastern_boundary_even;
	delete[] offset_to_eastern_boundary_odd;
	delete[] offset_to_northern_boundary_even;
	delete[] offset_to_northern_boundary_odd;
	delete[] offset_to_western_boundary_even;
	delete[] offset_to_western_boundary_odd;
	delete[] offset_to_southern_boundary_even;
	delete[] offset_to_southern_boundary_odd;
	delete[] inverse_dir;
	delete[] mirrored_dir_x;
	delete[] mirrored_dir_y;
	delete[] lattice_vec_x;
	delete[] lattice_vec_y;

	node_state_cpu                   = NULL;
	node_state_tmp_cpu               = NULL;
    node_state_out_cpu               = NULL;
	cell_type_cpu                    = NULL;
	cell_density_cpu                 = NULL;
	mean_density_cpu                 = NULL;
	cell_momentum_cpu                = NULL;
	mean_momentum_cpu                = NULL;
	offset_to_neighbor_even          = NULL;
	offset_to_neighbor_odd           = NULL;
	offset_to_eastern_boundary_even  = NULL;
	offset_to_eastern_boundary_odd   = NULL;
	offset_to_northern_boundary_even = NULL;
	offset_to_northern_boundary_odd  = NULL;
	offset_to_western_boundary_even  = NULL;
	offset_to_western_boundary_odd   = NULL;
	offset_to_southern_boundary_even = NULL;
	offset_to_southern_boundary_odd  = NULL;
	inverse_dir                      = NULL;
	mirrored_dir_x                   = NULL;
	mirrored_dir_y                   = NULL;
	lattice_vec_x                    = NULL;
	lattice_vec_y                    = NULL;
}

// Sets (proper) parallelization parameters.
void OMP_Lattice::setup_parallel()
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
std::vector<Real> OMP_Lattice::get_mean_velocity() {

    std::vector<Real> mean_velocity(dim, 0.0);

    Real sum_x_vel = 0.0;
    Real sum_y_vel = 0.0;

    unsigned int counter = 0;

    // Sum up all (fluid) cell x and y velocity components.
#pragma omp parallel for reduction(+: sum_x_vel, sum_y_vel)
    for (unsigned int n = 0; n < n_cells; ++n) {

        if (cell_type_cpu[n] == 0) {

        	counter++;

            Real cell_density = cell_density_cpu[n];

			if (cell_density > 1.0e-06) {

				sum_x_vel += cell_momentum_cpu[n          ] / cell_density;
				sum_y_vel += cell_momentum_cpu[n + n_cells] / cell_density;
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
