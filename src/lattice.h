/*
 * Lattice.h
 *
 *  Created on: Dec 9, 2015
 *      Author: Kerstin Vater
 * Description: This class defines a lattice gas cellular automaton in two
 *              dimensions.
 */

#ifndef LATTICE_H_
#define LATTICE_H_

// C++ include files
#include <vector>
#include <assert.h>

// PNGwriter include files
#include "pngwriter/pngwriter.h"

// User-defined include files
#include "utils.h"

// Namespaces
using namespace std;

// Define real number precision to use.
typedef float real;

class Lattice {

private:

protected:

    // Dimension of the problem.
    const int dim = 2;

    // Number of cells in x direction.
    unsigned int n_x;

    // Number of cells in y direction.
    unsigned int n_y;

    // Total number of cells.
    unsigned int n_cells;

    // Number of lattice directions defining the model under usage.
    int n_dir;

    // Total number of nodes in the lattice.
    unsigned int n_nodes;

    // Number of particles in the lattice.
    unsigned int n_particles;

    // Coarse graining radius, i.e. the number of neighbor cells
    // in one direction taken into account for averaging purposes.
    int coarse_graining_radius;

    // Test case (pipe flow, box, Karman vortex street, single collision).
    string test_case;

    // Reynolds number.
    real Re;

    // Scaled Mach number.
    real Ma_s;

    // Calculate the dimensions of the lattice from the specified Reynolds number.
    //
    // Density (is set to 1.0).
    const real rho = 1.0;

    // Speed of sound (is set to 1.0);
    const real c = 1.0;

    // Mean occupation number.
    real d;

    // Viscosity.
    real nu = 1.0 / 12.0 * 1.0 / (d * pow((1.0 - d), 3.0)) - 1.0 / 8.0;

    // Galilean breaking factor.
    real g;

    // Scaled viscosity.
    real nu_s;

    // Scaled sound speed.
    real c_s;

    // Velocity.
    real u;

    // Direction of the body force.
    char bf_dir;

    // Number of particles to revert in the context of body force
    // in order to compensate boundary layer shear force.
    int equilibrium_forcing;

    // Map which defines the type of the cells.
    //
    // 0 - fluid cell
    // 1 - solid cell, reflecting, bounce back
    // 2 - solid cell, reflecting, bounce forward
    char* cell_type_cpu;

    // One-dimensional arrays of integers which contains the states
    // of the nodes, i.e. the occupation numbers, of the cellular automaton
    // in the following sense:
    //
    // [DIR_1_CELL_1|DIR_1_CELL_2|DIR_1_CELL_3|...|DIR_2_CELL_1|DIR_2_CELL_2|...]
    //
    // Array on the CPU.
    char* node_state_cpu;

    // Density values (0th momentum) related to the single cells (non-averaged).
    real* cell_density_cpu;

    // Coarse grained density values (averaged over neighbor cells).
    real* mean_density_cpu;

    // Vector valued quantities are stored in one-dimensional arrays in the
    // following sense:
    //
    // [X_COMP_CELL_1|X_COMP_CELL_2|X_COMP_CELL_3|...|Y_COMP_CELL_1|Y_COMP_CELL_2|...]

    // Momentum vectors (1st momentum) related to the single cells (non-averaged).
    real* cell_momentum_cpu;

    // Coarse grained momentum vectors (averaged over neighbor cells).
    real* mean_momentum_cpu;

public:

    // Creates a lattice gas cellular automaton object of the specified properties.
    Lattice(const string test_case,
            const real Re, const real Ma_s,
            const int n_dir,
            const int coarse_graining_radius);

    // Deletes the lattice gas cellular automaton object.
    virtual ~Lattice();

    // Applies the specified cell type to all cells in the domain.
    void apply_cell_type_all(const int cell_type);

    // Applies the specified cell type to cells located on the eastern boundary
    // of the rectangular domain.
    void apply_boundary_cell_type_east(const int cell_type);

    // Applies the specified cell type to cells located on the northern boundary
    // of the rectangular domain.
    void apply_boundary_cell_type_north(const int cell_type);

    // Applies the specified cell type to cells located on the western boundary
    // of the rectangular domain.
    void apply_boundary_cell_type_west(const int cell_type);

    // Applies the specified cell type to cells located on the southern boundary
    // of the rectangular domain.
    void apply_boundary_cell_type_south(const int cell_type);

    // Applies boundary conditions for a peridoc domain, i.e. no boundaries over
    // the whole rectangular domain.
    void apply_bc_periodic();

    // Applies boundary conditions for reflecting boundaries at all edges
    // of the rectangular domain according to the specified bounce type
    // (back or forward).
    void apply_bc_reflecting(const string bounce_type);

    // Applies boundary conditions for a pipe flow, i.e. reflecting boundaries
    // at the upper and the lower edge of the rectangular domain.
    void apply_bc_pipe();

    // Applies boundary conditions for a Karman vortex street, i.e. a pipe flow
    // with a cylinder.
    void apply_bc_karman_vortex_street();

    // Initializes the lattice gas automaton with zero states.
    void init_zero();

    // Initializes the lattice gas automaton with some random distributed particles.
    void init_random();

    // Initializes the lattice gas automaton with single particles at defined nodes.
    void init_single(const vector<int> occupied_nodes);

    // Initializes the lattice gas automaton with two colliding particles.
    void init_single_collision();

    // Initializes the lattice gas automaton with some random distributed particles
    // in the center area of the domain.
    void init_diffusion();

    // Initializes the lattice gas automaton with some random distributed particles
    // in one quarter of a rectangular tank.
    void init_sloshing();

    // Returns the number of particles in the lattice.
    unsigned int get_n_particles();

    // Prints the lattice to the screen.
    void print();

    // Writes results of current time step to file.
    void write_results(const unsigned int step, const string format);

    // Prints information about the lattice object to screen.
    void print_info();

    // Computes the number of particles to revert in the context of body force
    // in order to compensate boundary layer shear force.
    int get_equilibrium_forcing();

    // Computes the number of particles to revert in the context of body force
    // in order to accelerate the flow.
    int get_initial_forcing();

    // Sets (proper) parallelization parameters.
    virtual void setup_parallel() = 0;

    // Calls the CUDA kernel which performs the collision and propagation step
    // on the lattice gas automaton.
    virtual void collide_and_propagate(unsigned int step) = 0;

    // Computes the mean velocity of the lattice.
    virtual vector<real> get_mean_velocity() = 0;

    // Calls the CUDA kernel which applies a body force in the specified
    // direction (x or y) and with the specified intensity to the particles.
    // E.g., if the intensity is equal 100, every 100th particle
    // changes it's direction, if feasible.
    virtual void apply_body_force(const int forcing) = 0;

    // Calls the CUDA kernels which compute quantities of interest as a
    // post-processing procedure.
    virtual void post_process() = 0;

    // Copies all data arrays from the host (CPU) to the device (GPU).
    virtual void copy_data_to_device();

    // Copies all data arrays from the device (GPU) back to the host (CPU).
    virtual void copy_data_from_device();

    // Get functions.
    real get_u()   { return u; }
    int  get_n_x() { return n_x; }
    int  get_n_y() { return n_y; }
};

#endif /* LATTICE_H_ */
