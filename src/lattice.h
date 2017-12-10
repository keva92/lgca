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

#ifndef LATTICE_H_
#define LATTICE_H_

#include "lgca_common.h"

#include "utils.h"

#include "pngwriter/pngwriter.h"

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
    Real Re;

    // Scaled Mach number.
    Real Ma_s;

    // Calculate the dimensions of the lattice from the specified Reynolds number.
    //
    // Density (is set to 1.0).
    const Real rho = 1.0;

    // Speed of sound (is set to 1.0);
    const Real c = 1.0;

    // Mean occupation number.
    Real d;

    // Viscosity.
    Real nu = 1.0 / 12.0 * 1.0 / (d * pow((1.0 - d), 3.0)) - 1.0 / 8.0;

    // Galilean breaking factor.
    Real g;

    // Scaled viscosity.
    Real nu_s;

    // Scaled sound speed.
    Real c_s;

    // Velocity.
    Real u;

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

    // One-dimensional arrays of integers which contains the states of the nodes, i.e. the
    // occupation numbers of the cellular automaton in the following sense:
    //
    // [DIR_1_CELL_1|DIR_1_CELL_2|DIR_1_CELL_3|...|DIR_2_CELL_1|DIR_2_CELL_2|...]
    //
    // Array on the CPU
    char* node_state_cpu;

    // Temporary buffer for post-processing and visualization
    char* node_state_out_cpu;

    // Density values (0th momentum) related to the single cells (non-averaged).
    Real* cell_density_cpu;

    // Coarse grained density values (averaged over neighbor cells).
    Real* mean_density_cpu;

    // Vector valued quantities are stored in one-dimensional arrays in the
    // following sense:
    //
    // [X_COMP_CELL_1|X_COMP_CELL_2|X_COMP_CELL_3|...|Y_COMP_CELL_1|Y_COMP_CELL_2|...]

    // Momentum vectors (1st momentum) related to the single cells (non-averaged).
    Real* cell_momentum_cpu;

    // Coarse grained momentum vectors (averaged over neighbor cells).
    Real* mean_momentum_cpu;

public:

    // Creates a lattice gas cellular automaton object of the specified properties.
    Lattice(const string test_case,
            const Real Re, const Real Ma_s,
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
    void init_single(const std::vector<int> occupied_nodes);

    // Initializes the lattice gas automaton with two colliding particles.
    void init_single_collision();

    // Initializes the lattice gas automaton with some random distributed particles
    // in the center area of the domain.
    void init_diffusion();

    // Returns the number of particles in the lattice.
    unsigned int get_n_particles();

    // Prints the lattice to the screen.
    void print();

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
    virtual std::vector<Real> get_mean_velocity() = 0;

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

    virtual void copy_data_to_output_buffer();

    // Get functions.
    Real get_u()   { return u; }
    int  get_n_x() { return n_x; }
    int  get_n_y() { return n_y; }

          Real* cell_density()       { assert(cell_density_cpu); return cell_density_cpu; }
    const Real* cell_density() const { assert(cell_density_cpu); return cell_density_cpu; }

          Real* mean_density()       { assert(mean_density_cpu); return mean_density_cpu; }
    const Real* mean_density() const { assert(mean_density_cpu); return mean_density_cpu; }

          Real* cell_momentum()       { assert(cell_momentum_cpu); return cell_momentum_cpu; }
    const Real* cell_momentum() const { assert(cell_momentum_cpu); return cell_momentum_cpu; }

          Real* mean_momentum()       { assert(mean_momentum_cpu); return mean_momentum_cpu; }
    const Real* mean_momentum() const { assert(mean_momentum_cpu); return mean_momentum_cpu; }

    Real cell_density(const int x, const int y) { assert(cell_density_cpu); return cell_density_cpu[y * n_x + x]; }
    Real mean_density(const int x, const int y) { assert(mean_density_cpu); return mean_density_cpu[y * n_x + x]; }
};

#endif /* LATTICE_H_ */
