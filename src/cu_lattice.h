/*
 * cu_lattice.h
 *
 *  Created on: Apr 4, 2016
 *      Author: Kerstin Vater
 * Description: This class defines a lattice gas cellular automaton in two
 *              dimensions parallelized by means of Nvidia CUDA.
 */

#ifndef CU_LATTICE_H_
#define CU_LATTICE_H_

#include "lgca_common.h"

#include "cuda_utils.cuh"
#include "lattice.h"

#include <curand.h>
#include <curand_kernel.h>

class CUDA_Lattice: public Lattice {

private:

    // Device to use for the simulation.
    int device;

    // Grid size in x, y and z direction.
    int grid_size_x;
    int grid_size_y;
    int grid_size_z;

    // Block size in x, y and z direction.
    int block_size_x;
    int block_size_y;
    int block_size_z;

    // Map which defines the type of the cells.
    //
    // 0 - fluid cell
    // 1 - solid cell, reflecting, bounce back
    // 2 - solid cell, reflecting, bounce forward
    char* cell_type_gpu;

    // One-dimensional arrays of integers which contains the states
    // of the nodes, i.e. the occupation numbers, of the cellular automaton
    // in the following sense:
    //
    // [DIR_1_CELL_1|DIR_1_CELL_2|DIR_1_CELL_3|...|DIR_2_CELL_1|DIR_2_CELL_2|...]
    //
    // Array on the GPU.
    char* node_state_gpu;

    // Auxiliary array on the GPU.
    char* node_state_tmp_gpu;

    // Density values (0th momentum) related to the single cells (non-averaged).
    Real* cell_density_gpu;

    // Coarse grained density values (averaged over neighbor cells).
    Real* mean_density_gpu;

    // Vector valued quantities are stored in one-dimensional arrays in the
    // following sense:
    //
    // [X_COMP_CELL_1|X_COMP_CELL_2|X_COMP_CELL_3|...|Y_COMP_CELL_1|Y_COMP_CELL_2|...]

    // Momentum vectors (1st momentum) related to the single cells (non-averaged).
    Real* cell_momentum_gpu;

    // Coarse grained momentum vectors (averaged over neighbor cells).
    Real* mean_momentum_gpu;

    // Allocates the memory for the arrays on the host (CPU) and device (GPU).
    void allocate_memory();

    // Frees the memory for the arrays on the host (CPU) and device (GPU).
    void free_memory();

    // Sets (proper) grid and block size for the GPU computation.
    void set_grid_and_block_size(int max_block_size);

public:

	CUDA_Lattice(const string test_case,
                 const Real Re, const Real Ma_s,
                 const int n_dir,
                 const int coarse_graining_radius,
                 const int device);

	virtual ~CUDA_Lattice();

    // Copies all data arrays from the host (CPU) to the device (GPU).
    void copy_data_to_device();

    // Copies all data arrays from the device (GPU) back to the host (CPU).
    void copy_data_from_device();

    // Sets (proper) parallelization parameters.
    void setup_parallel();

    // Calls the CUDA kernel which performs the collision and propagation step
    // on the lattice gas automaton.
    void collide_and_propagate(unsigned int step);

    // Computes the mean velocity of the lattice.
    std::vector<Real> get_mean_velocity();

    // Calls the CUDA kernel which applies a body force in the specified
    // direction (x or y) and with the specified intensity to the particles.
    // E.g., if the intensity is equal 100, every 100th particle
    // changes it's direction, if feasible.
    void apply_body_force(const int forcing);

    // Calls the CUDA kernels which computes quantities of interest as a
    // post-processing procedure.
    void post_process();
};

#endif /* CU_LATTICE_H_ */
