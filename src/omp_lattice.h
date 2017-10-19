/*
 * omp_lattice.h
 *
 *  Created on: Apr 4, 2016
 *      Author: Kerstin Vater
 * Description: This class defines a lattice gas cellular automaton in two
 *              dimensions parallelized by means of openMP.
 */

#ifndef OMP_LATTICE_H_
#define OMP_LATTICE_H_

#include "lattice.h"

class OMP_Lattice: public Lattice {

private:

    // Auxiliary array on the CPU.
    char* node_state_tmp_cpu;

	// Memory offset to neighbor cells in the different directions for the
	// propagation step.
	// Note that for the FHP model there is a difference in the offsets depending
	// on weather the cell is located in a row with even or odd index.
	int* offset_to_neighbor_even;
	int* offset_to_neighbor_odd;

	// Memory offset to related cells of the opposite boundary in the different
	// directions in case of periodic boundaries.
	int* offset_to_eastern_boundary_even;
	int* offset_to_eastern_boundary_odd;
	int* offset_to_northern_boundary_even;
	int* offset_to_northern_boundary_odd;
	int* offset_to_western_boundary_even;
	int* offset_to_western_boundary_odd;
	int* offset_to_southern_boundary_even;
	int* offset_to_southern_boundary_odd;

	// Inverse direction indices for each lattice direction.
	char* inverse_dir;

	// Mirrored direction indices for each lattice direction with respect
	// to the x and y axis.
	char* mirrored_dir_x;
	char* mirrored_dir_y;

	// Lattice vector components in the different directions.
	real* lattice_vec_x;
	real* lattice_vec_y;

	// Computes cell quantities of interest as a post-processing procedure.
	void cell_post_process();

	// Computes coarse grained quantities of interest as a post-processing procedure.
	void mean_post_process();

    // Allocates the memory for the arrays on the host (CPU) and device (GPU).
    void allocate_memory();

    // Frees the memory for the arrays on the host (CPU) and device (GPU).
    void free_memory();

public:

	// Creates an openMP parallelized lattice gas cellular automaton object
	// of the specified properties.
	OMP_Lattice(const string test_case,
                const real Re, const real Ma_s,
                const int n_dir,
                const int coarse_graining_radius);

	virtual ~OMP_Lattice();

    // Sets (proper) parallelization parameters.
    void setup_parallel();

    // Performs the collision and propagation step on the lattice gas automaton.
    void collide_and_propagate(unsigned int step);

    // Computes the mean velocity of the lattice.
    vector<real> get_mean_velocity();

    // Applies a body force in the specified direction (x or y) and with the
    // specified intensity to the particles. E.g., if the intensity is equal 100,
    // every 100th particle changes it's direction, if feasible.
    void apply_body_force(const int forcing);

    // Computes quantities of interest as a post-processing procedure.
    void post_process();
};

#endif /* OMP_LATTICE_H_ */
