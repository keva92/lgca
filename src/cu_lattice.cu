///*
// * cu_lattice.cu
// *
// *  Created on: Apr 4, 2016
// *      Author: Kerstin Vater
// * Description: This class defines a lattice gas cellular automaton in two
// *              dimensions parallelized by means of Nvidia CUDA.
// */

//#include "cu_lattice.h"

//__device__ Real cu_random(int seed, int offset) {

//    // CUDA's random number library uses curandState_t to keep
//    // track of the seed value we will store a random state for
//    // every thread.
//    curandState_t state;

//    // We have to initialize the state.
//    curand_init(seed,    // The seed controls the sequence of random values that are produced.
//                0,       // The sequence number is only important with multiple cores.
//                offset,  // The offset is how much extra we advance in the sequence for each call, can be 0.
//                &state);

//    return curand_uniform(&state);
//}

//__device__ bool cu_random_bool(int seed, int offset) {

//    return (cu_random(seed, offset) > 0.5);
//}

//// CUDA kernel which performs the collision and propagation step
//// on the lattice gas automaton.
//template <int n_dir>
//__global__ void collide_and_propagate_kernel(const int   n_x,
//                                             const int   n_y,
//                                             const int   n_cells,
//                                                   char* cell_type_gpu,
//                                                   char* node_state_gpu,
//                                                   char* node_state_tmp_gpu,
//                                             unsigned int step) {

//#ifdef DEBUG
//            // Check weather the domain dimensions are valid for the FHP model.
//            if (n_y % 2 != 0 && n_dir == 6) {

//                printf("ERROR in collide_and_propagate_kernel(): "
//                       "Invalid domain dimension in y direction.\n");
//                abort();
//            }
//#endif

//    // Each thread is working on one cell of the lattice.

//    // Thread is active per default.
//    bool active = true;

//    // Calculate the position of the cell in x direction (column index).
//    int pos_x = blockIdx.x * blockDim.x + threadIdx.x;

//    // Calculate the position of the cell in y direction (row index).
//    int pos_y = blockIdx.y;

//    // Check weather the thread is working on a valid cell.
//    if (pos_x >= n_x) {

//        active = false;
//    }

//    // Start calculation only for activated threads working on valid cells.
//    if (active) {

//        // Memory offset to neighbor cells in the different directions for the
//        // propagation step.
//        // Note that for the FHP model there is a difference in the offsets depending
//        // on weather the cell is located in a row with even or odd index.
//        int offset_to_neighbor[n_dir];

//        // Memory offset to related cells of the opposite boundary in the different
//        // directions in case of periodic boundaries.
//        int offset_to_eastern_boundary[n_dir];
//        int offset_to_northern_boundary[n_dir];
//        int offset_to_western_boundary[n_dir];
//        int offset_to_southern_boundary[n_dir];

//        // Inverse direction indices for each lattice direction.
//        char inverse_dir[n_dir];

//        // Mirrored direction indices for each lattice direction with respect
//        // to the x and y axis.
//        char mirrored_dir_x[n_dir];
//        char mirrored_dir_y[n_dir];

//        // Set the components of the lattice vectors for the different directions.
//        //
//        // Set the model based values according to the number of lattice directions.
//        switch (n_dir) {

//            // HPP model.
//            case 4:
//            {
//                offset_to_neighbor[0] = 1;
//                offset_to_neighbor[1] = n_x;
//                offset_to_neighbor[2] = -1;
//                offset_to_neighbor[3] = -n_x;

//                offset_to_eastern_boundary[0] = 0;
//                offset_to_eastern_boundary[1] = 0;
//                offset_to_eastern_boundary[2] = n_x;
//                offset_to_eastern_boundary[3] = 0;

//                offset_to_northern_boundary[0] = 0;
//                offset_to_northern_boundary[1] = 0;
//                offset_to_northern_boundary[2] = 0;
//                offset_to_northern_boundary[3] = n_x * n_y;

//                offset_to_western_boundary[0] = -n_x;
//                offset_to_western_boundary[1] = 0;
//                offset_to_western_boundary[2] = 0;
//                offset_to_western_boundary[3] = 0;

//                offset_to_southern_boundary[0] = 0;
//                offset_to_southern_boundary[1] = -n_x * n_y;
//                offset_to_southern_boundary[2] = 0;
//                offset_to_southern_boundary[3] = 0;

//                inverse_dir[0] = 2;
//                inverse_dir[1] = 3;
//                inverse_dir[2] = 0;
//                inverse_dir[3] = 1;

//                mirrored_dir_x[0] = 0;
//                mirrored_dir_x[1] = 3;
//                mirrored_dir_x[2] = 2;
//                mirrored_dir_x[3] = 1;

//                mirrored_dir_y[0] = 2;
//                mirrored_dir_y[1] = 1;
//                mirrored_dir_y[2] = 0;
//                mirrored_dir_y[3] = 3;

//                break;
//            }

//            // FHP model.
//            case 6:
//            {
//                // Define the memory offsets in the different directions for
//                // cells in rows with even and odd indices.
//                //
//                // The cell is located in a row with even index value.
//                if (pos_y % 2 == 0) {

//                    offset_to_neighbor[0] = 1;
//                    offset_to_neighbor[1] = n_x;
//                    offset_to_neighbor[2] = n_x - 1;
//                    offset_to_neighbor[3] = -1;
//                    offset_to_neighbor[4] = -n_x - 1;
//                    offset_to_neighbor[5] = -n_x;

//                    offset_to_eastern_boundary[0] = 0;
//                    offset_to_eastern_boundary[1] = 0;
//                    offset_to_eastern_boundary[2] = n_x;
//                    offset_to_eastern_boundary[3] = n_x;
//                    offset_to_eastern_boundary[4] = n_x;
//                    offset_to_eastern_boundary[5] = 0;

//                    offset_to_northern_boundary[0] = 0;
//                    offset_to_northern_boundary[1] = 0;
//                    offset_to_northern_boundary[2] = 0;
//                    offset_to_northern_boundary[3] = 0;
//                    offset_to_northern_boundary[4] = n_x * n_y;
//                    offset_to_northern_boundary[5] = n_x * n_y;

//                    offset_to_western_boundary[0] = -n_x;
//                    offset_to_western_boundary[1] = 0;
//                    offset_to_western_boundary[2] = 0;
//                    offset_to_western_boundary[3] = 0;
//                    offset_to_western_boundary[4] = 0;
//                    offset_to_western_boundary[5] = 0;

//                    offset_to_southern_boundary[0] = 0;
//                    offset_to_southern_boundary[1] = -n_x * n_y;
//                    offset_to_southern_boundary[2] = -n_x * n_y + 1;
//                    offset_to_southern_boundary[3] = 0;
//                    offset_to_southern_boundary[4] = 0;
//                    offset_to_southern_boundary[5] = 0;

//                // The cell is located in a row with odd index value.
//                } else if (pos_y % 2 != 0) {

//                    offset_to_neighbor[0] = 1;
//                    offset_to_neighbor[1] = n_x + 1;
//                    offset_to_neighbor[2] = n_x;
//                    offset_to_neighbor[3] = -1;
//                    offset_to_neighbor[4] = -n_x;
//                    offset_to_neighbor[5] = -n_x + 1;

//                    offset_to_eastern_boundary[0] = 0;
//                    offset_to_eastern_boundary[1] = 0;
//                    offset_to_eastern_boundary[2] = 0;
//                    offset_to_eastern_boundary[3] = n_x;
//                    offset_to_eastern_boundary[4] = 0;
//                    offset_to_eastern_boundary[5] = 0;

//                    offset_to_northern_boundary[0] = 0;
//                    offset_to_northern_boundary[1] = 0;
//                    offset_to_northern_boundary[2] = 0;
//                    offset_to_northern_boundary[3] = 0;
//                    offset_to_northern_boundary[4] = n_x * n_y;
//                    offset_to_northern_boundary[5] = n_x * n_y;

//                    offset_to_western_boundary[0] = -n_x;
//                    offset_to_western_boundary[1] = -n_x;
//                    offset_to_western_boundary[2] = 0;
//                    offset_to_western_boundary[3] = 0;
//                    offset_to_western_boundary[4] = 0;
//                    offset_to_western_boundary[5] = -n_x;

//                    offset_to_southern_boundary[0] = 0;
//                    offset_to_southern_boundary[1] = -n_x * n_y;
//                    offset_to_southern_boundary[2] = -n_x * n_y;
//                    offset_to_southern_boundary[3] = 0;
//                    offset_to_southern_boundary[4] = 0;
//                    offset_to_southern_boundary[5] = 0;
//                }

//                inverse_dir[0] = 3;
//                inverse_dir[1] = 4;
//                inverse_dir[2] = 5;
//                inverse_dir[3] = 0;
//                inverse_dir[4] = 1;
//                inverse_dir[5] = 2;

//                mirrored_dir_x[0] = 0;
//                mirrored_dir_x[1] = 5;
//                mirrored_dir_x[2] = 4;
//                mirrored_dir_x[3] = 3;
//                mirrored_dir_x[4] = 2;
//                mirrored_dir_x[5] = 1;

//                mirrored_dir_y[0] = 3;
//                mirrored_dir_y[1] = 2;
//                mirrored_dir_y[2] = 1;
//                mirrored_dir_y[3] = 0;
//                mirrored_dir_y[4] = 5;
//                mirrored_dir_y[5] = 4;

//                break;
//            }

//#ifdef DEBUG
//            default:
//            {
//                printf("ERROR in collide_and_propagate_kernel(): Invalid number of directions %d!\n", n_dir);
//                abort();
//                break;
//            }
//#endif

//        }

//        // Get index of the cell to work on.
//        int cell = n_x * blockIdx.y + pos_x;

//        // Get the type of the cell, i.e. fluid or solid.
//        // This has to be taken into account during the collision step, where
//        // cells behave different according to their type.
//        char cell_type = cell_type_gpu[cell];

//        // Check weather the cell is located on boundaries.
//        bool on_eastern_boundary  = ((blockIdx.x == (gridDim.x - 1)) && (pos_x == (n_x - 1)));
//        bool on_northern_boundary = (blockIdx.y == (gridDim.y - 1));
//        bool on_western_boundary  = ((threadIdx.x == 0) && (blockIdx.x == 0));
//        bool on_southern_boundary = (blockIdx.y == 0);

//        // Define an array for the global indices of the nodes in the cell.
//        int node_idx[n_dir];

//        // Define an array for the states of the nodes in the cell.
//        char node_state[n_dir];

//        // Execute collision step.
//        //
//        // The thread working on the cell has to know about the states of the
//        // nodes within the cell, therefore looping over all directions and
//        // look it up.
//#pragma unroll
//        for (int dir = 0; dir < n_dir; ++dir) {

//            node_idx[dir] = cell + dir * n_cells;
//            node_state[dir] = node_state_gpu[node_idx[dir]];
//        }

//        // TODO: Create a random boolean value for the collision step.
//        // bool rand_bool = cu_random_bool(seed, cell);
//        bool rand_bool =       ((pos_x % 2) == (pos_y % 2))
//                         - 1 * ((pos_x % 2) == (pos_y % 2)) * (step % 2)
//                         + 1 * ((pos_x % 2) != (pos_y % 2)) * (step % 2);

//        // Create a temporary array to copy the node states into.
//        char node_state_tmp[n_dir];

//        // Copy the actual states of the nodes to the temporary array.
//#pragma unroll
//        for (int dir = 0; dir < n_dir; ++dir) {

//            node_state_tmp[dir] = node_state[dir];
//        }

//        switch (cell_type) {

//            // The cell working on is a fluid cell ("normal" collision).
//            case 0:
//            {
//                // Using the the HPP model.
//                if (n_dir == 4) {

////                    // Collision case 1.
////                    if ((node_state[0] == 0) &&
////                        (node_state[1] == 1) &&
////                        (node_state[2] == 0) &&
////                        (node_state[3] == 1)) {
////
////                        node_state_tmp[0] = 1;
////                        node_state_tmp[1] = 0;
////                        node_state_tmp[2] = 1;
////                        node_state_tmp[3] = 0;
////
////                        break;
////                    }
////
////                    // Collision case 2.
////                    if ((node_state[0] == 1) &&
////                        (node_state[1] == 0) &&
////                        (node_state[2] == 1) &&
////                        (node_state[3] == 0)) {
////
////                        node_state_tmp[0] = 0;
////                        node_state_tmp[1] = 1;
////                        node_state_tmp[2] = 0;
////                        node_state_tmp[3] = 1;
////
////                        break;
////                    }

//                    node_state_tmp[0] = node_state[0]
//                            - (node_state[0] * node_state[2] * (1 - node_state[1]) * (1 - node_state[3]))
//                            + (node_state[1] * node_state[3] * (1 - node_state[0]) * (1 - node_state[2]));

//                    node_state_tmp[1] = node_state[1]
//                            - (node_state[1] * node_state[3] * (1 - node_state[0]) * (1 - node_state[2]))
//                            + (node_state[0] * node_state[2] * (1 - node_state[1]) * (1 - node_state[3]));

//                    node_state_tmp[2] = node_state[2]
//                            - (node_state[0] * node_state[2] * (1 - node_state[1]) * (1 - node_state[3]))
//                            + (node_state[1] * node_state[3] * (1 - node_state[0]) * (1 - node_state[2]));

//                    node_state_tmp[3] = node_state[3]
//                            - (node_state[1] * node_state[3] * (1 - node_state[0]) * (1 - node_state[2]))
//                            + (node_state[0] * node_state[2] * (1 - node_state[1]) * (1 - node_state[3]));

//                // Collision cases of the FHP model.
//                } else if (n_dir == 6) {

//                    // Collision case a1.
//                    if ((node_state[0] == 1) &&
//                        (node_state[1] == 0) &&
//                        (node_state[2] == 0) &&
//                        (node_state[3] == 1) &&
//                        (node_state[4] == 0) &&
//                        (node_state[5] == 0)) {

//                        node_state_tmp[0] = 0;
//                        node_state_tmp[1] = rand_bool;
//                        node_state_tmp[2] = 1 - node_state_tmp[1];
//                        node_state_tmp[3] = 0;
//                        node_state_tmp[4] = node_state_tmp[1];
//                        node_state_tmp[5] = node_state_tmp[2];

//                        break;
//                    }

//                    // Collision case a2.
//                    if ((node_state[0] == 0) &&
//                        (node_state[1] == 1) &&
//                        (node_state[2] == 0) &&
//                        (node_state[3] == 0) &&
//                        (node_state[4] == 1) &&
//                        (node_state[5] == 0)) {

//                        node_state_tmp[0] = rand_bool;
//                        node_state_tmp[1] = 0;
//                        node_state_tmp[2] = 1 - node_state_tmp[0];
//                        node_state_tmp[3] = node_state_tmp[0];
//                        node_state_tmp[4] = 0;
//                        node_state_tmp[5] = node_state_tmp[2];

//                        break;
//                    }

//                    // Collision case a3.
//                    if ((node_state[0] == 0) &&
//                        (node_state[1] == 0) &&
//                        (node_state[2] == 1) &&
//                        (node_state[3] == 0) &&
//                        (node_state[4] == 0) &&
//                        (node_state[5] == 1)) {

//                        node_state_tmp[0] = rand_bool;
//                        node_state_tmp[1] = 1 - node_state_tmp[0];
//                        node_state_tmp[2] = 0;
//                        node_state_tmp[3] = node_state_tmp[0];
//                        node_state_tmp[4] = node_state_tmp[1];
//                        node_state_tmp[5] = 0;

//                        break;
//                    }

//                    // Collision case b1.
//                    if ((node_state[0] == 0) &&
//                        (node_state[1] == 1) &&
//                        (node_state[2] == 0) &&
//                        (node_state[3] == 1) &&
//                        (node_state[4] == 0) &&
//                        (node_state[5] == 1)) {

//                        node_state_tmp[0] = 1;
//                        node_state_tmp[1] = 0;
//                        node_state_tmp[2] = 1;
//                        node_state_tmp[3] = 0;
//                        node_state_tmp[4] = 1;
//                        node_state_tmp[5] = 0;

//                        break;
//                    }

//                    // Collision case b2.
//                    if ((node_state[0] == 1) &&
//                        (node_state[1] == 0) &&
//                        (node_state[2] == 1) &&
//                        (node_state[3] == 0) &&
//                        (node_state[4] == 1) &&
//                        (node_state[5] == 0)) {

//                        node_state_tmp[0] = 0;
//                        node_state_tmp[1] = 1;
//                        node_state_tmp[2] = 0;
//                        node_state_tmp[3] = 1;
//                        node_state_tmp[4] = 0;
//                        node_state_tmp[5] = 1;

//                        break;
//                    }

////                    node_state_tmp[0] = node_state[0]
////                            - (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5]))
////                            + (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
////                            + (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4])) * (1 - rand_bool)
////                            - (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]))
////                            + (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]));
////
////                    node_state_tmp[1] = node_state[1]
////                            - (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5]))
////                            + (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
////                            + (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4])) * (1 - rand_bool)
////                            - (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]))
////                            + (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]));
////
////                    node_state_tmp[2] = node_state[2]
////                            - (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4]))
////                            + (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
////                            + (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5])) * (1 - rand_bool)
////                            - (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]))
////                            + (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]));
////
////                    node_state_tmp[3] = node_state[3]
////                            - (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5]))
////                            + (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
////                            + (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4])) * (1 - rand_bool)
////                            - (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]))
////                            + (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]));
////
////                    node_state_tmp[4] = node_state[4]
////                            - (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5]))
////                            + (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
////                            + (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4])) * (1 - rand_bool)
////                            - (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]))
////                            + (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]));
////
////                    node_state_tmp[5] = node_state[5]
////                            - (node_state[2] * node_state[5] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[1]) * (1 - node_state[4]))
////                            + (node_state[0] * node_state[3] * (1 - node_state[1]) * (1 - node_state[4]) * (1 - node_state[2]) * (1 - node_state[5])) * rand_bool
////                            + (node_state[1] * node_state[4] * (1 - node_state[0]) * (1 - node_state[3]) * (1 - node_state[2]) * (1 - node_state[5])) * (1 - rand_bool)
////                            - (node_state[1] * node_state[3] * node_state[5] * (1 - node_state[0]) * (1 - node_state[2]) * (1 - node_state[4]))
////                            + (node_state[0] * node_state[2] * node_state[4] * (1 - node_state[1]) * (1 - node_state[3]) * (1 - node_state[5]));
//                }

//#ifdef DEBUG
//                else {

//                    printf("ERROR in collide_and_propagate_kernel(): "
//                           "Invalid number of directions %d.\n", n_dir);
//                }
//#endif

//                break;
//            }

//            // The cell working on is a solid cell of bounce back type.
//            case 1:
//            {
//                // Loop over all directions.
//                // #pragma unroll
//                for (int dir = 0; dir < n_dir; ++dir) {

//                    // Exchange the states of the nodes with the the states of
//                    // the inverse directions.
//                    node_state_tmp[dir] = node_state[inverse_dir[dir]];
//                }

//                break;
//            }

//            // TODO: The cell working on is a solid cell of bounce forward type.
//            case 2:
//            {
//                // Loop over all directions.
//#pragma unroll
//                for (int dir = 0; dir < n_dir; ++dir) {

//                    if (on_northern_boundary || on_southern_boundary) {

//                        // Exchange the states of the nodes with the the states of
//                        // the mirrored directions along the x axis.
//                        node_state_tmp[dir] = node_state[mirrored_dir_x[dir]];
//                    }

//                    if (on_eastern_boundary || on_western_boundary) {

//                        // Exchange the states of the nodes with the the states of
//                        // the mirrored directions along the y axis.
//                        node_state_tmp[dir] = node_state[mirrored_dir_y[dir]];
//                    }
//                }

//                break;
//            }

//#ifdef DEBUG
//            // Invalid cell type.
//            default:
//            {
//                printf("ERROR in collide_and_propagate_kernel(): Invalid cell type %d.\n", cell_type);
//                break;
//            }
//#endif

//        }

//        // Execute propagation step.
//        //
//        // Loop over all directions.
//        // #pragma unroll
//        for (int dir = 0; dir < n_dir; dir++)
//        {
//            // Reset the memory offset.
//            int offset = 0;

//            // Construct the correct memory offset.
//            //
//            // Apply a default offset value.
//            offset += offset_to_neighbor[dir];

//            // Correct the offset in the current direction if the cell is
//            // located on boundaries.
//            if (on_eastern_boundary) {

//                offset += offset_to_western_boundary[dir];
//            }

//            if (on_northern_boundary) {

//                offset += offset_to_southern_boundary[dir];
//            }

//            if (on_western_boundary) {

//                offset += offset_to_eastern_boundary[dir];
//            }

//            if (on_southern_boundary) {

//                offset += offset_to_northern_boundary[dir];
//            }

//            // Push the states of the cell to its "neighbor" cells in the
//            // different directions.
//            node_state_tmp_gpu[node_idx[dir] + offset] = node_state_tmp[dir];
//        }
//    }
//}

//// Applies a body force in the specified direction (x or y) to the particles.
//template <int n_dir, char bf_dir>
//__global__ void apply_body_force_kernel(int   forcing,
//                                        int   n_x,
//                                        int   n_cells,
//                                        char* cell_type_gpu,
//                                        char* node_state_gpu,
//                                        char* node_state_tmp_gpu,
//                                        int   seed) {

//    // Each thread is looking for one particle to revert.

//    // Lattice vector components in the different directions.
//    Real lattice_vec_x[n_dir];
//    Real lattice_vec_y[n_dir];

//    // Mirrored direction indices for each lattice direction with respect
//    // to the x and y axis.
//    char mirrored_dir_x[n_dir];
//    char mirrored_dir_y[n_dir];

//    // Set the components of the lattice vectors for the different directions.
//    //
//    // Loop over all directions.
//    for (int dir = 0; dir < n_dir; ++dir) {

//        lattice_vec_x[dir] = cos(2.0 * M_PI / ((Real) n_dir) * ((Real) dir));
//        lattice_vec_y[dir] = sin(2.0 * M_PI / ((Real) n_dir) * ((Real) dir));
//    }

//    // Set the model based values according to the number of lattice directions.
//    switch (n_dir) {

//        // HPP model.
//        case 4:
//        {
//            mirrored_dir_x[0] = 0;
//            mirrored_dir_x[1] = 3;
//            mirrored_dir_x[2] = 2;
//            mirrored_dir_x[3] = 1;

//            mirrored_dir_y[0] = 2;
//            mirrored_dir_y[1] = 1;
//            mirrored_dir_y[2] = 0;
//            mirrored_dir_y[3] = 3;

//            break;
//        }

//        // FHP model.
//        case 6:
//        {
//            mirrored_dir_x[0] = 0;
//            mirrored_dir_x[1] = 5;
//            mirrored_dir_x[2] = 4;
//            mirrored_dir_x[3] = 3;
//            mirrored_dir_x[4] = 2;
//            mirrored_dir_x[5] = 1;

//            mirrored_dir_y[0] = 3;
//            mirrored_dir_y[1] = 2;
//            mirrored_dir_y[2] = 1;
//            mirrored_dir_y[3] = 0;
//            mirrored_dir_y[4] = 5;
//            mirrored_dir_y[5] = 4;

//            break;
//        }

//#ifdef DEBUG
//        default:
//        {
//            printf("ERROR in Lattice(): Invalid number of directions %d!\n", n_dir);
//            abort();
//            break;
//        }
//#endif

//    }

//    // Set a maximum number of iterations to find particles which can be reverted.
//    const int it_max = 2 * n_cells;

//    // Set the number of iterations to zero.
//    int it = 0;

//    // Number of particles which have been reverted.
//    int reverted_particles = 0;

//    // Loop over all cells.
//    do
//    {
//        // TODO: Get the index of a random cell.
//        int cell = (int)truncf(cu_random(seed, threadIdx.x) * n_cells);
//        printf("cell = %d", cell);
//        it++;

//        // Get the type of the cell, i.e. fluid or solid.
//        // Note that body forces are applied to fluid cells only.
//        char cell_type = cell_type_gpu[cell];

//        // Check weather the cell working on is a fluid cell.
//        if (cell_type == 0) {

//            // Define an array for the global indices of the nodes in the cell.
//            int node_idx[n_dir];

//            // Define an array for the states of the nodes in the cell.
//            char node_state[n_dir];

//            // The thread working on the cell has to know about the states of the
//            // nodes within the cell, therefore looping over all directions and
//            // look it up.
//        #pragma unroll
//            for (int dir = 0; dir < n_dir; ++dir) {

//                node_idx[dir] = cell + dir * n_cells;
//                node_state[dir] = node_state_gpu[node_idx[dir]];
//            }

//            // Create a temporary array to copy the node states into.
//            char node_state_tmp[n_dir];

//            // Copy the current states of the nodes to the temporary array.
//        #pragma unroll
//            for (int dir = 0; dir < n_dir; ++dir) {

//                node_state_tmp[dir] = node_state[dir];
//            }

//            if (n_dir == 4) {

//                if (bf_dir == 'x' && (node_state[0] == 0) && (node_state[2] == 1)) {

//                    node_state_tmp[0] = 1;
//                    node_state_tmp[2] = 0;

//                    reverted_particles++;

//                } else if (bf_dir == 'y' && (node_state[1] == 1) && (node_state[3] == 0)) {

//                    node_state_tmp[1] = 0;
//                    node_state_tmp[3] = 1;

//                    reverted_particles++;
//                }
//            }

//            else if (n_dir == 6) {

//                if (bf_dir == 'x' && (node_state[0] == 0) && (node_state[3] == 1)) {

//                    node_state_tmp[0] = 1;
//                    node_state_tmp[3] = 0;

//                    reverted_particles++;

//                } else if (bf_dir == 'y') {

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
//                }
//            }

//    //            // Loop over all directions.
//    //#pragma unroll
//    //            for (int dir = 0; dir < n_dir; ++dir) {
//    //
//    //                // Body force acting in x direction.
//    //                if (bf_dir == 'x') {
//    //
//    //					// TODO: Exchange the states of the nodes with the the states of
//    //					//       the mirrored directions along the y axis if feasible.
//    //					if ((fabs(lattice_vec_x[dir] - 1.0) < 1.0e-06) &&
//    //						(node_state[dir] < node_state[mirrored_dir_y[dir]])) {
//    //
//    //						node_state_tmp[dir                ] = node_state[mirrored_dir_y[dir]];
//    //						node_state_tmp[mirrored_dir_y[dir]] = node_state[dir                ];
//    //					}
//    //                }
//    //
//    //                // Body force acting in y direction.
//    //                else if (bf_dir == 'y') {
//    //
//    //					// TODO: Exchange the states of the nodes with the the states of
//    //					//       the mirrored directions along the x axis if feasible.
//    //					if ((lattice_vec_y[dir] < 1.0e-06) &&
//    //						(node_state[dir] < node_state[mirrored_dir_x[dir]])) {
//    //
//    //						node_state_tmp[dir                ] = node_state[mirrored_dir_x[dir]];
//    //						node_state_tmp[mirrored_dir_x[dir]] = node_state[dir                ];
//    //					}
//    //                }
//    //
//    //#ifdef DEBUG
//    //                // Invalid body force direction.
//    //                else {
//    //
//    //                    printf("ERROR in apply_body_force(): "
//    //                           "Invalid body force direction %c.\n", bf_dir);
//    //                }
//    //#endif
//    //            }

//            // Write the new node states back to the data array.
//            //
//            // Loop over all directions.
//        #pragma unroll
//            for(int dir = 0; dir < n_dir; dir++)
//            {
//                node_state_gpu[node_idx[dir]] = node_state_tmp[dir];
//            }

//        } /* IF cell_type */

//    } while ((reverted_particles < 1) && (it < it_max));
//}

//// Computes the mean velocity of the lattice.
//__global__ void get_mean_velocity_kernel(int *g_idata, int *g_odata) {

//    extern __shared__ int sdata[];

//    // Each thread loads one element from global to shared memory.
//    unsigned int tid = threadIdx.x;
//    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
//    sdata[tid] = g_idata[i];
//    __syncthreads();

//    // Do reduction in shared memory.
//    for (unsigned int s = 1; s < blockDim.x; s *= 2) {

//        if (tid % (2*s) == 0) {

//            sdata[tid] += sdata[tid + s];
//        }
//        __syncthreads();
//    }

//    // Write result for this block to global memory.
//    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
//}

//// Computes cell quantities of interest as a post-processing procedure.
//template <int n_dir>
//__global__ void cell_post_process_kernel(const int     n_x,
//                                         const int     n_cells,
//                                         const char*   node_state_gpu,
//                                         Real*         cell_density_gpu,
//                                         Real*         cell_momentum_gpu) {

//    // Each thread is working on one cell of the lattice.

//    // Thread is active per default.
//    bool active = true;

//    // Calculate the position of the cell in x direction.
//    int pos_x = blockIdx.x * blockDim.x + threadIdx.x;

//    // Check weather the thread is working on a valid cell.
//    if (pos_x >= n_x) {

//        active = false;
//    }

//    // Start calculation only for activated threads working on valid cells.
//    if (active) {

//        // Lattice vector components in the different directions.
//        Real lattice_vec_x[n_dir];
//        Real lattice_vec_y[n_dir];

//        // Set the components of the lattice vectors for the different directions.
//        //
//        // Loop over all directions.
//        for (int dir = 0; dir < n_dir; ++dir) {

//            lattice_vec_x[dir] = cos(2.0 * M_PI / ((Real) n_dir) * ((Real) dir));
//            lattice_vec_y[dir] = sin(2.0 * M_PI / ((Real) n_dir) * ((Real) dir));
//        }

//        // Get index of the cell to work on.
//        int cell = n_x * blockIdx.y + pos_x;

//        // Define an array for the global indices of the nodes in a cell.
//        int node_idx[n_dir];

//        // Define an array for the states of the nodes in a cell.
//        char node_state[n_dir];

//        // The thread working on the cell has to know about the states of the
//        // nodes within the cell, therefore looping over all directions and
//        // look it up.
//#pragma unroll
//        for (int dir = 0; dir < n_dir; ++dir) {

//            node_idx[dir] = cell + dir * n_cells;
//            node_state[dir] = node_state_gpu[node_idx[dir]];
//        }

//        // Initialize the cell quantities to be computed.
//        int  cell_density    = 0;
//        Real cell_momentum_x = 0.0;
//        Real cell_momentum_y = 0.0;

//        // Loop over all nodes in the cell.
//#pragma unroll
//        for (int dir = 0; dir < n_dir; ++dir) {

//            // Sum up the node states.
//            cell_density += node_state[dir];

//            // Sum up the node states multiplied by the lattice vector component
//            // for the current direction.
//            cell_momentum_x += node_state[dir] * lattice_vec_x[dir];
//            cell_momentum_y += node_state[dir] * lattice_vec_y[dir];
//        }

//        // Write the computed cell quantities to the related data arrays.
//        cell_density_gpu [cell          ] = (Real) cell_density;
//        cell_momentum_gpu[cell          ] =        cell_momentum_x;
//        cell_momentum_gpu[cell + n_cells] =        cell_momentum_y;
//    }
//}

//// Computes coarse grained quantities of interest as a post-processing procedure.
//__global__ void mean_post_process_kernel(const int   n_x,
//                                         const int   n_cells,
//                                         const int   coarse_graining_radius,
//                                         const Real* cell_density_gpu,
//                                         const Real* cell_momentum_gpu,
//                                               Real* mean_density_gpu,
//                                               Real* mean_momentum_gpu) {

//    // Each thread is working on one cell of the lattice.

//    // Thread is active per default.
//    bool active = true;

//    // Calculate the position of the cell in x direction.
//    int pos_x = blockIdx.x * blockDim.x + threadIdx.x;

//    // Check weather the thread is working on a valid cell.
//    if (pos_x >= n_x) {

//        active = false;
//    }

//    // Start calculation only for activated threads working on valid cells.
//    if (active) {

//        // Get index of the cell to work on.
//        int cell = n_x * blockIdx.y + pos_x;

//        // Initialize the coarse grained quantities to be computed.
//        Real mean_density    = 0.0;
//        Real mean_momentum_x = 0.0;
//        Real mean_momentum_y = 0.0;

//        // Initialize the number of actual existing coarse graining neighbor cells.
//        int n_exist_neighbors = 0;

//        // The thread working on the cell has to know the cell quantities of the
//        // coarse graining neighbor cells, therefore looping over all neighbor
//        // cells and look it up.
//#pragma unroll
//        for (int y = -coarse_graining_radius; y <= coarse_graining_radius; ++y) {

//            for (int x = -coarse_graining_radius; x <= coarse_graining_radius; ++x) {

//                // Get the index of the coarse graining neighbor cell.
//                int neighbor_idx = cell + y * n_x + x;

//                // Get the position of the coarse graining neighbor cell in x direction.
//                int pos_x_neighbor = neighbor_idx % n_x;

//                // Check weather the coarse graining neighbor cell is valid.
//                if ((neighbor_idx >= 0) &&
//                    (neighbor_idx < n_cells) &&
//                    (abs(pos_x_neighbor - pos_x) <= coarse_graining_radius)) {

//                    // Increase the number of existing coarse graining neighbor cells.
//                    n_exist_neighbors++;

//                    mean_density    += cell_density_gpu [neighbor_idx          ];
//                    mean_momentum_x += cell_momentum_gpu[neighbor_idx          ];
//                    mean_momentum_y += cell_momentum_gpu[neighbor_idx + n_cells];
//                }
//            }
//        }

//        // Write the computed coarse grained quantities to the related data arrays.
//        mean_density_gpu [cell          ] = mean_density    / ((Real) n_exist_neighbors);
//        mean_momentum_gpu[cell          ] = mean_momentum_x / ((Real) n_exist_neighbors);
//        mean_momentum_gpu[cell + n_cells] = mean_momentum_y / ((Real) n_exist_neighbors);
//    }
//}

//// Creates a CUDA parallelized lattice gas cellular automaton object
//// of the specified properties.
//CUDA_Lattice::CUDA_Lattice(const string test_case,
//                           const Real Re, const Real Ma_s,
//                           const int n_dir,
//                           const int coarse_graining_radius,
//                           const int device = 0)

//                           : Lattice(test_case, Re, Ma_s, n_dir, coarse_graining_radius) {

//    // Set the device to use for the simulation.
//    int n_devices;
//    cudaGetDeviceCount(&n_devices);
//    assert((device < n_devices) && (device >= 0));
//    this->device = device;
//    cudaSetDevice(device);

//    // Allocate the memory for the arrays on the host (CPU) and device (GPU).
//    allocate_memory();
//}

//// Deletes the CUDA parallelized lattice gas cellular automaton object.
//CUDA_Lattice::~CUDA_Lattice() {

//    free_memory();
//}

//// Sets (proper) grid and block sizes for the GPU computation.
//void CUDA_Lattice::set_grid_and_block_size(int max_block_size = 256) {

//    grid_size_x = 1;
//    grid_size_y = n_y;
//    grid_size_z = 1;

//    block_size_x = n_x;
//    block_size_y = 1;
//    block_size_z = 1;

//    if (n_x > max_block_size) {

//        grid_size_x  = (int) (ceil((Real)n_x / (Real)max_block_size) + 0.5);
//        block_size_x = max_block_size;
//        if (n_x % max_block_size != 0) {
//            printf("WARNING in Lattice::set_grid_and_block_size(): "
//                   "There are inactive threads in some blocks!\n");
//        }
//    }

//    cudaDeviceProp prop;
//    cudaGetDeviceProperties(&prop, device);

//    if ( grid_size_x                               <=   prop.maxGridSize[0] &&
//         grid_size_y                               <=   prop.maxGridSize[1] &&
//         grid_size_z                               <=   prop.maxGridSize[2] &&
//        block_size_x                               <= prop.maxThreadsDim[0] &&
//        block_size_y                               <= prop.maxThreadsDim[1] &&
//        block_size_z                               <= prop.maxThreadsDim[2] &&
//        block_size_x * block_size_y * block_size_z <= prop.maxThreadsPerBlock) {

//        printf("Kernel configuration parameters: %d x %d x %d Blocks \n",   grid_size_x , grid_size_y , grid_size_z );
//        printf("                                 %d x %d x %d Threads\n\n", block_size_x, block_size_y, block_size_z);

//    } else {

//        printf("ERROR in Lattice::set_grid_and_block_size():"
//               "Invalid grid and/or block dimensions. "
//               "Please check device properties.");
//        abort();
//    }
//}

//// Copies all data arrays from the device (GPU) back to the host (CPU).
//void CUDA_Lattice::copy_data_from_device() {

//    cu_verify(cudaMemcpy(node_state_cpu,    node_state_gpu,          n_nodes * sizeof(char), cudaMemcpyDeviceToHost));
//    cu_verify(cudaMemcpy(cell_density_cpu,  cell_density_gpu,        n_cells * sizeof(Real), cudaMemcpyDeviceToHost));
//    cu_verify(cudaMemcpy(mean_density_cpu,  mean_density_gpu,        n_cells * sizeof(Real), cudaMemcpyDeviceToHost));
//    cu_verify(cudaMemcpy(cell_momentum_cpu, cell_momentum_gpu, dim * n_cells * sizeof(Real), cudaMemcpyDeviceToHost));
//    cu_verify(cudaMemcpy(mean_momentum_cpu, mean_momentum_gpu, dim * n_cells * sizeof(Real), cudaMemcpyDeviceToHost));
//}

//// Copies all data arrays from the host (CPU) to the device (GPU).
//void CUDA_Lattice::copy_data_to_device() {

//    cu_verify(cudaMemcpy(cell_type_gpu,  cell_type_cpu,  n_cells * sizeof(char), cudaMemcpyHostToDevice));
//    cu_verify(cudaMemcpy(node_state_gpu, node_state_cpu, n_nodes * sizeof(char), cudaMemcpyHostToDevice));
//}

//// Allocates the memory for the arrays on the host (CPU) and device (GPU).
//void CUDA_Lattice::allocate_memory() {

//    // Allocate host memory.
//    cu_verify(cudaMallocHost((void **) &node_state_cpu,          n_nodes * sizeof(char)));
//    cu_verify(cudaMallocHost((void **) &cell_type_cpu,           n_cells * sizeof(char)));
//    cu_verify(cudaMallocHost((void **) &cell_density_cpu,        n_cells * sizeof(Real)));
//    cu_verify(cudaMallocHost((void **) &mean_density_cpu,        n_cells * sizeof(Real)));
//    cu_verify(cudaMallocHost((void **) &cell_momentum_cpu, dim * n_cells * sizeof(Real)));
//    cu_verify(cudaMallocHost((void **) &mean_momentum_cpu, dim * n_cells * sizeof(Real)));

//    // Allocate device memory.
//    cu_verify(cudaMalloc((void **) &node_state_gpu,          n_nodes * sizeof(char)));
//    cu_verify(cudaMalloc((void **) &node_state_tmp_gpu,      n_nodes * sizeof(char)));
//    cu_verify(cudaMalloc((void **) &cell_type_gpu,           n_cells * sizeof(char)));
//    cu_verify(cudaMalloc((void **) &cell_density_gpu,        n_cells * sizeof(Real)));
//    cu_verify(cudaMalloc((void **) &mean_density_gpu,        n_cells * sizeof(Real)));
//    cu_verify(cudaMalloc((void **) &cell_momentum_gpu, dim * n_cells * sizeof(Real)));
//    cu_verify(cudaMalloc((void **) &mean_momentum_gpu, dim * n_cells * sizeof(Real)));
//}

//// Frees the memory for the arrays on the host (CPU) and device (GPU).
//void CUDA_Lattice::free_memory() {

//    // Free GPU memory.
//    cu_verify(cudaFree(node_state_gpu));
//    cu_verify(cudaFree(node_state_tmp_gpu));
//    cu_verify(cudaFree(cell_type_gpu));
//    cu_verify(cudaFree(cell_density_gpu));
//    cu_verify(cudaFree(mean_density_gpu));
//    cu_verify(cudaFree(cell_momentum_gpu));
//    cu_verify(cudaFree(mean_momentum_gpu));

//    // Free CPU memory.
//    cu_verify(cudaFreeHost(node_state_cpu));
//    cu_verify(cudaFreeHost(cell_type_cpu));
//    cu_verify(cudaFreeHost(cell_density_cpu));
//    cu_verify(cudaFreeHost(mean_density_cpu));
//    cu_verify(cudaFreeHost(cell_momentum_cpu));
//    cu_verify(cudaFreeHost(mean_momentum_cpu));
//}

//// Calls the CUDA kernel which performs the collision and propagation step
//// on the lattice gas automaton.
//void CUDA_Lattice::collide_and_propagate(unsigned int step) {

//    // Set the grid and block size.
//    dim3 grid_size (grid_size_x,  grid_size_y,  grid_size_z);
//    dim3 block_size(block_size_x, block_size_y, block_size_z);

//    // TODO: Set the seed for the random number generation on the device.
//    int seed = time(NULL);

//    // Call CUDA kernel.
//    switch (n_dir) {

//        // HPP model.
//        case 4:
//        {
//            cu_verify_void((collide_and_propagate_kernel<4>
//                    <<<grid_size, block_size>>>(n_x,
//                                                n_y,
//                                                n_cells,
//                                                cell_type_gpu,
//                                                node_state_gpu,
//                                                node_state_tmp_gpu,
//                                                step)));
//            break;
//        }

//        // FHP model.
//        case 6:
//        {
//            cu_verify_void((collide_and_propagate_kernel<6>
//                    <<<grid_size, block_size>>>(n_x,
//                                                n_y,
//                                                n_cells,
//                                                cell_type_gpu,
//                                                node_state_gpu,
//                                                node_state_tmp_gpu,
//                                                step)));
//            break;
//        }
//        default:
//        {
//            printf("ERROR in collide_and_propagate(): "
//                   "Invalid number of directions %d.\n", n_dir);
//            abort();

//            break;
//        }
//    }

//    // Wait for all threads to finish.
//    cudaDeviceSynchronize();

//    // Update the node states.
//    char* node_state_gpu_tmp = node_state_gpu;
//    node_state_gpu = node_state_tmp_gpu;
//    node_state_tmp_gpu = node_state_gpu_tmp;
//}

//// Calls the CUDA kernel which applies a body force in the specified
//// direction (x or y) to the particles.
//void CUDA_Lattice::apply_body_force(const int forcing) {

//    const int max_block_size = 256;

//    // Set the grid and block size.
//    int grid_size_x = forcing / max_block_size + 1;

//    int block_size_x = 1;
//    if (forcing > 0)   block_size_x = forcing;
//    if (forcing > 256) block_size_x = max_block_size;

//    dim3 grid_size (grid_size_x,  1, 1);
//    dim3 block_size(block_size_x, 1, 1);

//    // Get device properties.
//    cudaDeviceProp prop;
//    cudaGetDeviceProperties(&prop, device);

//    if ( grid_size.x                               <= prop.maxGridSize[0] &&
//         grid_size.y                               <= prop.maxGridSize[1] &&
//         grid_size.z                               <= prop.maxGridSize[2] &&
//        block_size.x                               <= prop.maxThreadsDim[0] &&
//        block_size.y                               <= prop.maxThreadsDim[1] &&
//        block_size.z                               <= prop.maxThreadsDim[2] &&
//        block_size.x * block_size.y * block_size.z <= prop.maxThreadsPerBlock) {

//        printf("Kernel configuration parameters: %d x %d x %d Blocks \n", grid_size.x , grid_size.y , grid_size.z );
//        printf("                                 %d x %d x %d Threads\n", block_size.x, block_size.y, block_size.z);

//    } else {

//        printf("ERROR in CUDA_Lattice::apply_body_force():"
//               "Invalid grid and/or block dimensions. "
//               "Please check device properties.");
//        abort();
//    }

//    // TODO: Set the seed for the random number generation on the device.
//    int seed = time(NULL);

//    // Call CUDA kernel.
//    switch (n_dir) {

//        // HPP model.
//        case 4:
//        {
//            switch (bf_dir) {

//                // Apply body force in x direction.
//                case 'x':
//                {
//                    cu_verify_void((apply_body_force_kernel<4, 'x'>
//                            <<<grid_size, block_size>>>(forcing,
//                                                        n_x,
//                                                        n_cells,
//                                                        cell_type_gpu,
//                                                        node_state_gpu,
//                                                        node_state_tmp_gpu,
//                                                        seed)));
//                    break;
//                }

//                // Apply body force in y direction.
//                case 'y':
//                {
//                    cu_verify_void((apply_body_force_kernel<4, 'y'>
//                            <<<grid_size, block_size>>>(forcing,
//                                                        n_x,
//                                                        n_cells,
//                                                        cell_type_gpu,
//                                                        node_state_gpu,
//                                                        node_state_tmp_gpu,
//                                                        seed)));
//                    break;
//                }

//                // Invalid body force direction.
//                default:
//                {
//                    printf("ERROR in apply_body_force(): "
//                           "Invalid body force direction %c.\n", bf_dir);
//                    abort();

//                    break;
//                }
//            }

//            break;
//        }

//        // FHP model.
//        case 6:
//        {
//            switch (bf_dir) {

//                // Apply body force in x direction.
//                case 'x':
//                {
//                    cu_verify_void((apply_body_force_kernel<6, 'x'>
//                            <<<grid_size, block_size>>>(forcing,
//                                                        n_x,
//                                                        n_cells,
//                                                        cell_type_gpu,
//                                                        node_state_gpu,
//                                                        node_state_tmp_gpu,
//                                                        seed)));
//                    break;
//                }

//                // Apply body force in y direction.
//                case 'y':
//                {
//                    cu_verify_void((apply_body_force_kernel<6, 'y'>
//                            <<<grid_size, block_size>>>(forcing,
//                                                        n_x,
//                                                        n_cells,
//                                                        cell_type_gpu,
//                                                        node_state_gpu,
//                                                        node_state_tmp_gpu,
//                                                        seed)));
//                    break;
//                }

//                // Invalid body force direction.
//                default:
//                {
//                    printf("ERROR in apply_body_force(): "
//                           "Invalid body force direction %c.\n", bf_dir);
//                    abort();

//                    break;
//                }
//            }

//            break;
//        }

//        // Invalid number of directions.
//        default:
//        {
//            printf("ERROR in apply_body_force(): "
//                   "Invalid number of directions %d.\n", n_dir);
//            abort();

//            break;
//        }
//    }


//    // Wait for all threads to finish.
//    cudaDeviceSynchronize();
//}

//// Call the CUDA kernel which computes quantities of interest as a
//// post-processing procedure.
//void CUDA_Lattice::post_process() {

//    // Set the grid and block size.
//    dim3 grid_size (grid_size_x,  grid_size_y,  grid_size_z);
//    dim3 block_size(block_size_x, block_size_y, block_size_z);

//    // Call CUDA kernel.
//    switch (n_dir) {

//        // HPP model.
//        case 4:
//        {
//            cu_verify_void((cell_post_process_kernel<4>
//                    <<<grid_size, block_size>>>(n_x,
//                                                n_cells,
//                                                node_state_gpu,
//                                                cell_density_gpu,
//                                                cell_momentum_gpu)));
//            break;
//        }

//        // FHP model.
//        case 6:
//        {
//            cu_verify_void((cell_post_process_kernel<6>
//                    <<<grid_size, block_size>>>(n_x,
//                                                n_cells,
//                                                node_state_gpu,
//                                                cell_density_gpu,
//                                                cell_momentum_gpu)));
//            break;
//        }

//        // Invalid number of directions.
//        default:
//        {
//            printf("ERROR in post_process(): "
//                   "Invalid number of directions %d.\n", n_dir);
//            abort();

//            break;
//        }
//    }

//    // Wait for all threads to finish.
//    cudaDeviceSynchronize();

//    cu_verify_void((mean_post_process_kernel
//            <<<grid_size, block_size>>>(n_x,
//                                        n_cells,
//                                        coarse_graining_radius,
//                                        cell_density_gpu,
//                                        cell_momentum_gpu,
//                                        mean_density_gpu,
//                                        mean_momentum_gpu)));

//    // Wait for all threads to finish.
//    cudaDeviceSynchronize();
//}

//// Sets (proper) parallelization parameters.
//void CUDA_Lattice::setup_parallel()
//{
//    // Sets (proper) grid and block size for the GPU computation.
//    set_grid_and_block_size(256);
//}

//// TODO: Computes the mean velocity of the lattice.
//std::vector<Real> CUDA_Lattice::get_mean_velocity()
//{
//    std::vector<Real> mean_velocity(dim, 0.0);

//    Real sum_x_vel = 0.0;
//    Real sum_y_vel = 0.0;

//    unsigned int counter = 0;

//    // Sum up all (fluid) cell x and y velocity components.
//#pragma omp parallel for reduction(+: sum_x_vel, sum_y_vel)
//    for (unsigned int n = 0; n < n_cells; ++n) {

//        if (cell_type_cpu[n] == 0) {

//            counter++;

//            Real cell_density = cell_density_cpu[n];

//            if (cell_density > 1.0e-06) {

//                sum_x_vel += cell_momentum_cpu[n          ] / cell_density;
//                sum_y_vel += cell_momentum_cpu[n + n_cells] / cell_density;
//            }

//#ifdef DEBUG

//            else if (fabs(cell_density) < 1.0e-06) {

//                // Do nothing.

//            } else if (cell_density < -1.0e-06) {

//                printf("ERROR in get_mean_velocity(): "
//                       "Negative cell density detected.");
//                abort();
//            }

//#endif

//        }
//    }

//    // Divide the summed up x and y components by the total number of fluid cells.
//    mean_velocity[0] = sum_x_vel / (Real) counter;
//    mean_velocity[1] = sum_y_vel / (Real) counter;

//    return mean_velocity;
//}

