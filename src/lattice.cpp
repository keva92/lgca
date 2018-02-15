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

#include "lattice.h"

#include <omp.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace lgca {

// Creates a lattice gas cellular automaton object of the specified properties.
template<Model model_>
Lattice<model_>::Lattice(const string test_case,
                         const Real Re, const Real Ma_s,
                         const int coarse_graining_radius) {

    // Set the test case (pipe flow,
	//					  box,
	//					  Karman vortex street,
	//                    single collision,
    //                    diffusion).
    this->m_test_case = test_case;

    // Set the Reynolds number.
    assert(Re > 1.0e-06);
    this->m_Re = Re;

    // Set the Mach number.
    assert(Ma_s > 1.0e-06);
    this->m_Ma_s = Ma_s;

    // Calculate the dimensions of the lattice from the specified Reynolds number.
    //
    // Define the mean occupation number.
    m_d = m_rho / NUM_DIR;

    // Compute the scaled viscosity.
    m_nu = 1.0 / 12.0 * 1.0 / (m_d * pow((1.0 - m_d), 3.0)) - 1.0 / 8.0;

    // Compute the Galilean breaking factor.
    m_g = SPATIAL_DIM / (SPATIAL_DIM + 2.0) * (1.0 - 2.0 * m_d) / (1.0 - m_d);

    // Compute the viscosity.
    m_nu_s = m_nu / m_g;

    // Get the scaled sound speed.
    m_c_s = m_c / sqrt(SPATIAL_DIM);

    // Get the velocity.
    m_u = this->m_Ma_s * m_c_s;

    // Get the number of cells in y direction;
    if (test_case == "pipe") {

        // Get the number of cells in y direction in case of a pipe flow.
        m_dim_y = (int)((Re * m_nu_s) / m_u + 0.5);

    } else if (test_case == "karman") {

        // Get the cylinder diameter in case of a Karman vortex street.
        Real diameter = (Re * m_nu_s) / m_u;
        m_dim_y = (int)(3.0 * diameter + 0.5);

    } else if (test_case == "collision") {

        // Set the number of cells in y direction to 8.
        m_dim_y = 8;

    } else if (test_case == "diffusion" ||
    		   test_case == "periodic"  ||
               test_case == "box") {

        m_dim_y = (int)Re;

    } else {

        printf("ERROR in Lattice::Lattice(): Invalid test case %s.\n", test_case.c_str());
        abort();
    }

    // Define the number of cells in x direction.
    if (test_case == "pipe"      ||
        test_case == "karman") {

        m_dim_x = 2 * m_dim_y;

    } else if (test_case == "collision" ||
               test_case == "box"       ||
    	       test_case == "diffusion" ||
    	       test_case == "periodic") {

        m_dim_x = m_dim_y;

	} else {

		printf("ERROR in Lattice::Lattice(): Invalid test case %s.\n", test_case.c_str());
		abort();
	}

    // Correct number of cells to fit the bitset container type
    // (automatically matches FHP requirement for even number of cells in y direction)
    const size_t cell_blocks_wide = (m_dim_x + Bitset::BITS_PER_BLOCK-1) / Bitset::BITS_PER_BLOCK + 2;
    const size_t cell_blocks_high = (m_dim_y +                        1) / Bitset::BITS_PER_BLOCK + 1;
    m_dim_x = cell_blocks_wide * Bitset::BITS_PER_BLOCK;
    m_dim_y = cell_blocks_high * Bitset::BITS_PER_BLOCK;

    // Set the body force direction according to the test case.
    if (test_case == "pipe"   ||
        test_case == "karman") {

        m_bf_dir = 'x';

    } else {

        m_bf_dir = 0;
    }

    // Set the number of cells in x direction.
    assert(m_dim_x > 0);
    this->m_dim_x = m_dim_x;

    // Set the number of cells in y direction.
    assert(m_dim_y > 0);
    if ((model_ == Model::FHP_I) || (model_ == Model::FHP_II) || (model_ == Model::FHP_III)) assert(m_dim_y % 2 == 0);
    this->m_dim_y = m_dim_y;

    m_num_cells = m_dim_x * m_dim_y;
    m_num_nodes = m_num_cells * NUM_DIR;

    m_num_particles = 0;

    assert(coarse_graining_radius > 0);
    this->m_coarse_graining_radius = coarse_graining_radius;

    this->m_coarse_dim_x     = m_dim_x / (2 * m_coarse_graining_radius); assert(m_dim_x % (2 * m_coarse_graining_radius) == 0);
    this->m_coarse_dim_y     = m_dim_y / (2 * m_coarse_graining_radius); assert(m_dim_y % (2 * m_coarse_graining_radius) == 0);
    this->m_num_coarse_cells = m_coarse_dim_x * m_coarse_dim_y;

    print_info();
}

// Deletes the lattice gas cellular automaton object.
template<Model model_>
Lattice<model_>::~Lattice() {

}

// Initializes the lattice gas automaton with zero states.
template<Model model_>
void Lattice<model_>::init_zero() {

//    memset(m_node_state_cpu, 0, m_num_nodes * sizeof(unsigned char));
}

// Prints the lattice to the screen.
template<Model model_>
void Lattice<model_>::print() {

    m_node_state_cpu.print();
}

// Returns the number of particles in the lattice.
template<Model model_>
unsigned long Lattice<model_>::get_n_particles() {

    // TODO Implement and use count() function in Bitset class

    size_t n_particles = 0;

    // Loop over all the nodes.
#pragma omp parallel for reduction(+:n_particles)
    for (size_t n = 0; n < m_num_nodes; ++n)
        n_particles += bool(m_node_state_cpu[n]);

    this->m_num_particles = n_particles;

    return n_particles;
}

// Initializes the lattice gas automaton with some random distributed particles.
template<Model model_>
void Lattice<model_>::init_random() {

    // Loop over all cells
#pragma omp parallel for
    for (size_t cell = 0; cell < m_num_cells; ++cell) {

        // Check weather the cell is a fluid cell
        if (m_cell_type_cpu[cell] == bool(CellType::FLUID)) {

            // Loop over all nodes in the fluid cell.
            for (int dir = 0; dir < NUM_DIR; ++dir) {

                // Set random states for the nodes in the fluid cell
                m_node_state_cpu[cell + dir * m_num_cells] =
                        bool(random_uniform() > (1.0 - (1.0 / NUM_DIR)));
            }
	    }
	}
}

// Applies boundary conditions for a peridoc domain, i.e. no boundaries over
// the whole rectangular domain.
template<Model model_>
void Lattice<model_>::apply_bc_periodic() {

    // Set the cell type to fluid cells
    CellType cell_type = CellType::FLUID;

    // Set all cell types to fluid cells
    apply_cell_type_all(cell_type);
}

// Applies boundary conditions for a pipe flow, i.e. reflecting boundaries
// at the upper and the lower edge of the rectangular domain.
template<Model model_>
void Lattice<model_>::apply_bc_pipe() {

    // Set the cell type to fluid cells
    CellType cell_type = CellType::FLUID;

    // Set all cell types to fluid cells
    apply_cell_type_all(cell_type);

    // Set the cell type to solid boundary cells of bounce back type
    cell_type = CellType::SOLID_NO_SLIP;

    // Apply the specified cell type to cells located on the boundaries of the rectangular domain
    apply_boundary_cell_type_north(cell_type);
    apply_boundary_cell_type_south(cell_type);
}

// Applies boundary conditions for a Karman vortex street, i.e. a pipe flow
// with a cylinder.
template<Model model_>
void Lattice<model_>::apply_bc_karman_vortex_street() {

    // Apply boundary conditions for a pipe flow, i.e. reflecting boundaries at the upper and the
    // lower edge of the rectangular domain
    apply_bc_pipe();

    // Define the position and size of the barrier
    int  center_x = m_dim_x / 6;
    int  center_y = m_dim_y / 2 + 1 / 10 * m_dim_y;
    Real diameter = m_dim_y / 3;

    // Loop over all cells
#pragma omp parallel for
    for (size_t cell = 0; cell < m_num_cells; ++cell) {

        // Get the position of the current cell
        int pos_x = cell % m_dim_x;
        int pos_y = cell / m_dim_x;

        Real dist = sqrt(pow((pos_x - center_x), 2.0) + pow((pos_y - center_y), 2.0));

        if (dist < (diameter / 2.0)) {

            // Set the cell type to solid cells of bounce back type
            m_cell_type_cpu[cell] = bool(CellType::SOLID_NO_SLIP);
        }
    }
}

// Initializes the lattice gas automaton with two colliding particles
template<Model model_>
void Lattice<model_>::init_single_collision() {

    std::vector<size_t> occupied_nodes;
    occupied_nodes.push_back((1 * m_dim_x + m_dim_x / 2) * NUM_DIR +                    0);
    occupied_nodes.push_back((3 * m_dim_x + m_dim_x / 2) * NUM_DIR + ModelDesc::INV_DIR[0]);

    init_single(occupied_nodes);
}

// Initializes the lattice gas automaton with single particles at defined nodes
template<Model model_>
void Lattice<model_>::init_single(const std::vector<size_t> occupied_nodes) {

    // Initialize the lattice with zeros
    init_zero();

    // Loop over the nodes to place a particle at
    for (int n = 0; n < occupied_nodes.size(); ++n) {

        m_node_state_cpu[occupied_nodes[n]] = bool(1);
    }
}

// Applies boundary conditions for reflecting boundaries at all edges of the rectangular domain
template<Model model_>
void Lattice<model_>::apply_bc_reflecting(const string bounce_type) {

    // Set the cell type to fluid cells
    CellType cell_type = CellType::FLUID;

    // Set all cell types to fluid cells
    apply_cell_type_all(cell_type);

    // Get the type of the boundary cells according to the bounce type
    if (bounce_type == "back") cell_type = CellType::SOLID_NO_SLIP;
    else {

        printf("ERROR in apply_bc_reflecting(): Invalid bounce type %s.\n", bounce_type.c_str());
        abort();
    }

    // Apply the specified cell type to cells located on the boundaries of the rectangular domain
    apply_boundary_cell_type_east (cell_type);
    apply_boundary_cell_type_north(cell_type);
    apply_boundary_cell_type_west (cell_type);
    apply_boundary_cell_type_south(cell_type);
}

// Applies the specified cell type to all cells in the domain
template<Model model_>
void Lattice<model_>::apply_cell_type_all(const CellType cell_type) {

    // Loop over all cells
#pragma omp parallel for
    for (size_t cell = 0; cell < m_num_cells; ++cell) {

        // Set the cell type to the specified type
        m_cell_type_cpu[cell] = bool(cell_type);
    }
}

// Applies the specified cell type to cells located on the eastern boundary of the rectangular domain
template<Model model_>
void Lattice<model_>::apply_boundary_cell_type_east(const CellType cell_type) {

    // Loop over the cells located at the eastern boundary of the rectangular domain
    for (size_t cell = m_dim_x - 1; cell < m_num_cells; cell += m_dim_x) {

        // Set the cell type to the specified type
        m_cell_type_cpu[cell] = bool(cell_type);
    }
}

// Applies the specified cell type to cells located on the northern boundary of the rectangular domain
template<Model model_>
void Lattice<model_>::apply_boundary_cell_type_north(const CellType cell_type) {

    // Loop over the cells located at the northern boundary of the rectangular domain
    for (size_t cell = m_num_cells - m_dim_x; cell < m_num_cells; ++cell) {

        // Set the cell type to the specified type
        m_cell_type_cpu[cell] = bool(cell_type);
    }
}

// Applies the specified cell type to cells located on the western boundary of the rectangular domain
template<Model model_>
void Lattice<model_>::apply_boundary_cell_type_west(const CellType cell_type) {

    // Loop over the cells located at the western boundary of the rectangular domain
    for (size_t cell = 0; cell < m_num_cells; cell += m_dim_x) {

        // Set the cell type to the specified type
        m_cell_type_cpu[cell] = bool(cell_type);
    }
}

// Applies the specified cell type to cells located on the southern boundary of the rectangular domain
template<Model model_>
void Lattice<model_>::apply_boundary_cell_type_south(const CellType cell_type) {

    // Loop over the cells located at the southern boundary of the rectangular domain
    for (size_t cell = 0; cell < m_dim_x; ++cell) {

        // Set the cell type to the specified type
        m_cell_type_cpu[cell] = bool(cell_type);
    }
}

// Prints information about the lattice object to screen.
template<Model model_>
void Lattice<model_>::print_info() {

    printf("Parameter for test case \"%s\":\n", m_test_case.c_str());
    printf("\n");
    printf("Reynolds number             Re   = %10.2f\n", m_Re);
    printf("Mach number                 Ma   = %10.2f\n", m_Ma_s);
    printf("Density                     rho  = %10.2f\n", m_rho);
    printf("Velocity                    u    = %10.2f\n", m_u);
    printf("Sound speed                 c    = %10.2f\n", m_c);
    printf("Scaled sound speed          c_s  = %10.2f\n", m_c_s);
    printf("Viscosity                   nu   = %10.2f\n", m_nu);
    printf("Scaled viscosity            nu_s = %10.2f\n", m_nu_s);
    printf("Galilean breaking factor    g    = %10.2f\n", m_g);
    printf("\n");
    printf("Number of cells in x direction: %d\n", m_dim_x);
    printf("Number of cells in y direction: %d\n", m_dim_y);
    printf("\n");
    printf("Number of coarse cells in x direction: %d\n", m_coarse_dim_x);
    printf("Number of coarse cells in y direction: %d\n", m_coarse_dim_y);
    printf("\n");
}

// Copies all data arrays from the host (CPU) to the device (GPU).
template<Model model_>
void Lattice<model_>::copy_data_to_device()
{
	// Only valid for CUDA parallelized lattice gas automatons.
	// Implemented in CUDA_Lattice.
}

// Copies all data arrays from the device (GPU) back to the host (CPU).
template<Model model_>
void Lattice<model_>::copy_data_from_device()
{
	// Only valid for CUDA parallelized lattice gas automatons.
	// Implemented in CUDA_Lattice.
}

template<Model model_>
void Lattice<model_>::copy_data_to_output_buffer()
{
    m_node_state_out_cpu.copy(m_node_state_cpu);
}

// Computes the number of particles to revert in the context of body force
// in order to accelerate the flow.
template<Model model_>
size_t Lattice<model_>::get_initial_forcing()
{
    // size_t equilibrium_forcing = get_equilibrium_forcing();

    return (size_t)(0.01 * m_num_cells);
}

// Computes the number of particles to revert in the context of body force.
// in order to compensate boundary layer shear force.
template<Model model_>
size_t Lattice<model_>::get_equilibrium_forcing()
{
    Real forcing = (8.0 * m_nu_s * m_Ma_s * m_c_s) / pow((Real)m_dim_y, 2.0);

    return ceil(0.5 * m_num_cells * forcing);
}

// Initializes the lattice gas automaton with some random distributed particles
// in the center area of the domain.
template<Model model_>
void Lattice<model_>::init_diffusion()
{
    // Define the position and size of the center area
    const int  center_x = m_dim_x / 2;
    const int  center_y = m_dim_y / 2;
    const Real diameter = m_dim_y / 4;

    Bitset::Block* __restrict__ write = m_node_state_cpu.ptr();

    // Loop over bunches of cells
    const size_t num_cell_blocks = ((m_num_cells - 1) / Bitset::BITS_PER_BLOCK) + 1;
    const size_t cell_block_dim_x = m_dim_x / Bitset::BITS_PER_BLOCK;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_cell_blocks), [&](const tbb::blocked_range<size_t>& r) {
    for (size_t cell_block = r.begin(); cell_block != r.end(); ++cell_block)
    {
        // Get the y position of the current cell block
        const size_t cell_block_pos_y = cell_block / cell_block_dim_x;

        Bitset outputs(Bitset::BITS_PER_BLOCK * NUM_DIR);

        // Loop over cells in cell block
#pragma unroll
        for (int local_cell = 0; local_cell < Bitset::BITS_PER_BLOCK; ++local_cell) {

            const size_t global_cell = cell_block * Bitset::BITS_PER_BLOCK + local_cell;

            // Get the x and y position of the current cell
            const int cell_pos_x = global_cell % m_dim_x;
            const int cell_pos_y = cell_block_pos_y;

            Real dist = sqrt(pow((cell_pos_x - center_x), 2.0) + pow((cell_pos_y - center_y), 2.0));

            // Check weather the cell is a fluid cell in the center area of the domain
            if (m_cell_type_cpu[global_cell] == bool(CellType::FLUID) && dist < (diameter / 2.0))
            {
                // Loop over all nodes in the fluid cell
                for (int dir = 0; dir < NUM_DIR; ++dir) {

                    // Set random states for the nodes in the fluid cell
                    outputs[local_cell + dir * Bitset::BITS_PER_BLOCK] =
                            bool(random_uniform() > (1.0 - (1.0 / NUM_DIR)));
                }
            }

        } // for local cell

#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir)
            write[cell_block * NUM_DIR + dir] = outputs(dir);

    }}); // for cell block
}

template<Model model_>
void Lattice<model_>::init_pipe() {

    Bitset::Block* __restrict__ write = m_node_state_cpu.ptr();

    // Loop over bunches of cells
    const size_t num_cell_blocks = ((m_num_cells - 1) / Bitset::BITS_PER_BLOCK) + 1;
    const size_t cell_block_dim_x = m_dim_x / Bitset::BITS_PER_BLOCK;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_cell_blocks), [&](const tbb::blocked_range<size_t>& r) {
    for (size_t cell_block = r.begin(); cell_block != r.end(); ++cell_block)
    {
        // Get the y position of the current cell block
        const size_t cell_block_pos_y = cell_block / cell_block_dim_x;

        const Real y = 1.0 * cell_block_pos_y / m_dim_y;
        const Real factor = 6.0 * (1.0 - y) * y; // TODO Check
        const Real u_y = 0.0;
        const Real u_x = factor * m_u;

        Bitset outputs(Bitset::BITS_PER_BLOCK * NUM_DIR);

        // Loop over cells in cell block
#pragma unroll
        for (int local_cell = 0; local_cell < Bitset::BITS_PER_BLOCK; ++local_cell) {

            const size_t global_cell = cell_block * Bitset::BITS_PER_BLOCK + local_cell;

            // Check weather the cell is a fluid cell
            if (m_cell_type_cpu[global_cell] == bool(CellType::FLUID)) {

                // Loop over all nodes in the fluid cell
#pragma unroll
                for (int dir = 0; dir < NUM_DIR; ++dir) {

                    // Calculate equilibrium distribution function for direction dir
                    Real N_eq = m_rho / NUM_DIR + m_rho * 2.0 / NUM_DIR * (ModelDesc::LATTICE_VEC_X[dir] * u_x
                                                                         + ModelDesc::LATTICE_VEC_Y[dir] * u_y);

                    // Set random states for the nodes in the fluid cell
                    outputs[local_cell + dir * Bitset::BITS_PER_BLOCK] =
                            bool(random_uniform() < N_eq);
                }
            }

        } // for local cell

#pragma unroll
        for (int dir = 0; dir < NUM_DIR; ++dir)
            write[cell_block * NUM_DIR + dir] = outputs(dir);

    }}); // for cell block
}

// Explicit instantiations
template class Lattice<Model::HPP>;
template class Lattice<Model::FHP_I>;
template class Lattice<Model::FHP_II>;
template class Lattice<Model::FHP_III>;

} // namespace lgca
