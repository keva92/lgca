/*
 * lattice.cpp
 *
 *  Created on: Dec 9, 2015
 *      Author: Kerstin Vater
 * Description: This class defines a lattice gas cellular automaton in two
 *              dimensions.
 *
 */

#include "lattice.h"

#include <omp.h>

// Creates a lattice gas cellular automaton object of the specified properties.
Lattice::Lattice(const string test_case,
                 const Real Re, const Real Ma_s,
                 const int n_dir,
                 const int coarse_graining_radius) {

    // Set the test case (pipe flow,
	//					  box,
	//					  Karman vortex street,
	//                    single collision,
	//                    diffusion,
	//                    sloshing,
	//                    hourglass).
    this->test_case = test_case;

    // Set the number of lattice directions defining the model under usage.
    assert((n_dir == 4) || (n_dir == 6));
    this->n_dir = n_dir;

    // Set the Reynolds number.
    assert(Re > 1.0e-06);
    this->Re = Re;

    // Set the Mach number.
    assert(Ma_s > 1.0e-06);
    this->Ma_s = Ma_s;

    // Calculate the dimensions of the lattice from the specified Reynolds number.
    //
    // Define the mean occupation number.
    d = rho / this->n_dir;

    // Compute the scaled viscosity.
    nu = 1.0 / 12.0 * 1.0 / (d * pow((1.0 - d), 3.0)) - 1.0 / 8.0;

    // Compute the Galilean breaking factor.
    g = dim / (dim + 2.0) * (1.0 - 2.0 * d) / (1.0 - d);

    // Compute the viscosity.
    nu_s = nu / g;

    // Get the scaled sound speed.
    c_s = c / sqrt(dim);

    // Get the velocity.
    u = this->Ma_s * c_s;

    // Get the number of cells in y direction;
    if (test_case == "pipe") {

        // Get the number of cells in y direction in case of a pipe flow.
        n_y = (int)((Re * nu_s) / u + 0.5);

    } else if (test_case == "karman") {

        // Get the cylinder diameter in case of a Karman vortex street.
        Real diameter = (Re * nu_s) / u;
        n_y = (int)(3.0 * diameter + 0.5);

    } else if (test_case == "collision") {

        // Set the number of cells in y direction to 5.
        n_y = 7;

    } else if (test_case == "diffusion" ||
    		   test_case == "periodic"  ||
    		   test_case == "box"       ||
    		   test_case == "sloshing"  ||
    		   test_case == "hourglass") {

    	n_y = (int)Re;

    } else {

        printf("ERROR in Lattice::Lattice(): Invalid test case %s.\n", test_case.c_str());
        abort();
    }

    // Correct the number of cells in y direction for use with the FHP mode.
    if (n_y % 2 != 0) n_y++;

    // Define the number of cells in x direction.
    if (test_case == "pipe"      ||
        test_case == "karman"    ||
        test_case == "sloshing"  ||
        test_case == "hourglass" ||
        test_case == "collision") {

    	n_x = 2 * n_y;

    } else if (test_case == "box"       ||
    	       test_case == "diffusion" ||
    	       test_case == "periodic") {

    	n_x = n_y;

	} else {

		printf("ERROR in Lattice::Lattice(): Invalid test case %s.\n", test_case.c_str());
		abort();
	}

    // Correct the number of cells in x direction for use with collision debug test case.
    if (n_x % 2 != 0) n_x++;
    if (test_case == "collision") n_x++;

    // Set the body force direction according to the test case.
    if (test_case == "pipe"   ||
        test_case == "karman" ||
        test_case == "hourglass") {

    	bf_dir = 'x';

    } else if (test_case == "sloshing") {

    	bf_dir = 'y';

    } else {

        bf_dir = 0;
    }

    // Set the number of cells in x direction.
    assert(n_x > 0);
    this->n_x = n_x;

    // Set the number of cells in y direction.
    assert(n_y > 0);
    if(n_dir == 6) assert(n_y % 2 == 0);
    this->n_y = n_y;

    // Set the total number of cells.
    n_cells = n_x * n_y;

    // Set the total number of nodes in the lattice.
    n_nodes = n_cells * n_dir;

    // Set the initial number of particles in the lattice.
    n_particles = 0;

    // Set the coarse graining radius.
    assert(coarse_graining_radius > 0);
    this->coarse_graining_radius = coarse_graining_radius;

    print_info();
}

// Deletes the lattice gas cellular automaton object.
Lattice::~Lattice() {

}

// Initializes the lattice gas automaton with zero states.
void Lattice::init_zero() {

    // Loop over all the nodes.
#pragma omp parallel for
    for (int n = 0; n < n_nodes; ++n) {

        // Initialize the lattice with zeros.
        node_state_cpu[n] = 0;
    }
}

// Prints the lattice to the screen.
void Lattice::print() {

    printf("lattice = [\n");

    // Loop over all the nodes.
    for (int n = 0; n < n_nodes; ++n) {

        printf("%d ", node_state_cpu[n]);

        if (((n + 1) % n_x == 0) && (n != 0)) {

            printf("\n");
        }

        if (((n + 1) % n_cells == 0) && (n != 0)) {

            printf("\n");
        }
    }

    printf("]\n");

    printf("\n");

    printf("cell_type = [\n");

    // Loop over all the cells.
    for (int n = 0; n < n_cells; ++n) {

        printf("%d ", cell_type_cpu[n]);

        if (((n + 1) % n_x == 0) && (n != 0)) {

            printf("\n");
        }
    }

    printf("]\n");
}

// Returns the number of particles in the lattice.
unsigned int Lattice::get_n_particles() {

    int n_particles = 0;

    // Loop over all the nodes.
#pragma omp parallel for reduction(+:n_particles)
    for (unsigned int n = 0; n < n_nodes; ++n) {

        n_particles += node_state_cpu[n];
    }

    this->n_particles = n_particles;

    return n_particles;
}

// Initializes the lattice gas automaton with some random distributed particles.
void Lattice::init_random() {

	// Loop over all cells.
#pragma omp parallel for
	for (int cell = 0; cell < n_cells; ++cell) {

	    // Check weather the cell is a fluid cell.
	    if (cell_type_cpu[cell] == 0) {

            // Loop over all nodes in the fluid cell.
            for (int dir = 0; dir < n_dir; ++dir) {

                // Set random states for the nodes in the fluid cell.
            	node_state_cpu[cell + dir * n_cells] =
            	        (random_uniform() > (1.0 - (1.0 / n_dir)));
            }
	    }
	}
}

// Applies boundary conditions for a peridoc domain, i.e. no boundaries over
// the whole rectangular domain.
void Lattice::apply_bc_periodic() {

    // Set the cell type to fluid cells.
    int cell_type = 0;

    // Set all cell types to fluid cells.
    apply_cell_type_all(cell_type);
}

// Applies boundary conditions for a pipe flow, i.e. reflecting boundaries
// at the upper and the lower edge of the rectangular domain.
void Lattice::apply_bc_pipe() {

    // Set the cell type to fluid cells.
    int cell_type = 0;

    // Set all cell types to fluid cells.
    apply_cell_type_all(cell_type);

    // Set the cell type to solid boundary cells of bounce back type.
    cell_type = 1;

    // Apply the specified cell type to cells located on the boundaries
    // of the rectangular domain.
    apply_boundary_cell_type_north(cell_type);
    apply_boundary_cell_type_south(cell_type);
}

// Applies boundary conditions for a Karman vortex street, i.e. a pipe flow
// with a cylinder.
void Lattice::apply_bc_karman_vortex_street() {

    // Apply boundary conditions for a pipe flow, i.e. reflecting boundaries
    // at the upper and the lower edge of the rectangular domain.
    apply_bc_pipe();

    // Define the position and size of the barrier.
    int  center_x = n_x / 6;
    int  center_y = n_y / 2 + 1 / 10 * n_y;
    Real diameter = n_y / 3;

    // Loop over all cells.
#pragma omp parallel for
    for (int cell = 0; cell < n_cells; ++cell) {

        // Get the position of the current cell.
        int pos_x = cell % n_x;
        int pos_y = cell / n_x;

        Real dist = sqrt(pow((pos_x - center_x), 2.0) + pow((pos_y - center_y), 2.0));

        if (dist < (diameter / 2.0)) {

            // Set the cell type to solid cells of bounce back type.
            cell_type_cpu[cell] = 1;
        }
    }
}

// Initializes the lattice gas automaton with two colliding particles.
void Lattice::init_single_collision() {

	// Get the inverse direction of the 0-th one.
	int inverse_dir;

	// HPP model.
	if (n_dir == 4) {

		inverse_dir = 2;

	// FHP model.
	} else if (n_dir == 6) {

		inverse_dir = 3;

	// Invalid number of directions.
	} else {

        printf("ERROR in init_single_collision(): Invalid number of directions %d!\n", n_dir);
        abort();
	}

    std::vector<int> occupied_nodes;
    occupied_nodes.push_back(n_x * n_y / 4 + 1);
    occupied_nodes.push_back(occupied_nodes[0] + (n_x - 9) + inverse_dir * n_cells);

    init_single(occupied_nodes);
}

// Writes results of current time step to file.
void Lattice::write_results(const unsigned int step, const string format) {

    // Set flags which data should be written to file.
    bool write_cell_density  = true;
    bool write_mean_density  = true;
    bool write_cell_momentum = false;
    bool write_mean_momentum = false;
    bool write_cell_velocity = false;
    bool write_mean_velocity = false;

	// Create a file.
	FILE* file;
	char name[100], end[100], s_iter[100];

	// Write vti file.
	if (format == "vti") {

	    // Create file name.
        sprintf(name, "../res/res_");
	    sprintf(s_iter, "%d", step);
	    strcpy(end, format.c_str());
	    strcat(name, s_iter);
	    strcat(name, ".");
	    strcat(name, end);

	    printf("Writing results to file %s...\n", name);

        // Open the result file.
        file = fopen(name, "w");
        assert(file);

        fprintf(file, "<VTKFile "
             "type=\"ImageData\" "
             "version=\"0.1\" "
             "byte_order=\"LittleEndian\" "
             "compressor=\"vtkZLibDataCompressor\">\n");

        fprintf(file, "  <ImageData "
             "WholeExtent=\"0 %d 0 %d 0 0\" "
             "Origin=\"0 0 0\" "
             "Spacing=\"1.0 1.0 1.0\">\n",
             (n_x - 1), (n_y - 1));

        fprintf(file, "  <Piece Extent="
             "\"0 %d 0 %d 0 0\">\n",
             (n_x - 1), (n_y - 1));

        fprintf(file, "    <PointData>\n");

        if (write_cell_density == true) {

            // Write cell density values from lattice to file.
            fprintf(file, "      <DataArray "
                 "type=\"Float32\" "
                 "Name=\"Cell Density\" "
                 "format=\"ascii\">\n", n_dir);

            // Loop over all cells.
            for (int cell = 0; cell < n_cells; ++cell) {

                fprintf(file, "      %f\n", cell_density_cpu[cell]);
            }

            fprintf(file, "      </DataArray>\n");
        }

        if (write_mean_density == true) {

            // Write mean density values from lattice to file.
            fprintf(file, "      <DataArray "
                 "type=\"Float32\" "
                 "Name=\"Mean Density\" "
                 "format=\"ascii\">\n", n_dir);

            // Loop over all cells.
            for (int cell = 0; cell < n_cells; ++cell) {

                if (mean_density_cpu[cell] > -1.0e-06) {

                    fprintf(file, "      %f\n", mean_density_cpu[cell]);
                }
#ifdef DEBUG
                else {

                    printf("ERROR in write_results(): "
                           "Negative mean density detected.");
                    abort();
                }
#endif
            }

            fprintf(file, "      </DataArray>\n");
        }

        if (write_cell_momentum == true) {

            // Write cell momentum vectors from lattice to file.
            fprintf(file, "      <DataArray "
                 "type=\"Float32\" "
                 "Name=\"Cell Momentum\" "
                 "NumberOfComponents=\"3\" "
                 "format=\"ascii\">\n");

            // Loop over all cells.
            for (int cell = 0; cell < n_cells; ++cell) {

                fprintf(file, "      ");

                fprintf(file, "%f ", cell_momentum_cpu[cell]);
                fprintf(file, "%f ", cell_momentum_cpu[cell + n_cells]);
                fprintf(file, "%f ", 0.0); // Print zero as third component.

                fprintf(file, "\n");
            }

            fprintf(file, "      </DataArray>\n");
        }

        if (write_mean_momentum == true) {

            // Write mean momentum vectors from lattice to file.
            fprintf(file, "      <DataArray "
                 "type=\"Float32\" "
                 "Name=\"Mean Momentum\" "
                 "NumberOfComponents=\"3\" "
                 "format=\"ascii\">\n");

            // Loop over all cells.
            for (int cell = 0; cell < n_cells; ++cell) {

                fprintf(file, "      ");

                fprintf(file, "%f ", mean_momentum_cpu[cell]);
                fprintf(file, "%f ", mean_momentum_cpu[cell + n_cells]);
                fprintf(file, "%f ", 0.0); // Print zero as third component.

                fprintf(file, "\n");
            }

            fprintf(file, "      </DataArray>\n");
        }

        if (write_cell_velocity == true) {

            // Write cell velocity vectors from lattice to file.
            fprintf(file, "      <DataArray "
                 "type=\"Float32\" "
                 "Name=\"Cell Velocity\" "
                 "NumberOfComponents=\"3\" "
                 "format=\"ascii\">\n");

            // Loop over all cells.
            for (int cell = 0; cell < n_cells; ++cell) {

                fprintf(file, "      ");

                if (cell_density_cpu[cell] > 1.0e-06) {

                    fprintf(file, "%f ", cell_momentum_cpu[cell          ] / cell_density_cpu[cell]);
                    fprintf(file, "%f ", cell_momentum_cpu[cell + n_cells] / cell_density_cpu[cell]);

                } else if (fabs(cell_density_cpu[cell]) < 1.0e-06) {

                    fprintf(file, "%f ", 0.0);
                    fprintf(file, "%f ", 0.0);
                }

        #ifdef DEBUG

                else {

                    printf("ERROR in write_results(): "
                           "Negative cell density detected.");
                    abort();
                }

        #endif

                fprintf(file, "%f ", 0.0); // Print zero as third component.

                fprintf(file, "\n");
            }

            fprintf(file, "      </DataArray>\n");
        }

        if (write_mean_velocity == true) {

            // Write mean velocity vectors from lattice to file.
            fprintf(file, "      <DataArray "
                 "type=\"Float32\" "
                 "Name=\"Mean Velocity\" "
                 "NumberOfComponents=\"3\" "
                 "format=\"ascii\">\n");

            // Loop over all cells.
            for (int cell = 0; cell < n_cells; ++cell) {

                fprintf(file, "      ");

                if (mean_density_cpu[cell] > 1.0e-06) {

                    fprintf(file, "%f ", mean_momentum_cpu[cell]           / mean_density_cpu[cell]);
                    fprintf(file, "%f ", mean_momentum_cpu[cell + n_cells] / mean_density_cpu[cell]);

                } else if (fabs(mean_density_cpu[cell]) < 1-0e-06) {

                    fprintf(file, "%f ", 0.0);
                    fprintf(file, "%f ", 0.0);

                }

        #ifdef DEBUG

                else {

                    printf("ERROR in write_results(): "
                           "Negative mean density detected.");
                    abort();
                }

        #endif

                fprintf(file, "%f ", 0.0); // Print zero as third component.

                fprintf(file, "\n");
            }

            fprintf(file, "      </DataArray>\n");
        }

        fprintf(file, "    </PointData>\n");

        fprintf(file, "    <CellData>\n");

        fprintf(file, "    </CellData>\n");

        fprintf(file, "  </Piece>\n");

        fprintf(file, "  </ImageData>\n");

        fprintf(file, "</VTKFile>\n");

        // Close the result file.
        fclose(file);

    // Write png file.
	} else if (format == "png") {

        if (write_mean_density == true) {

            // Create file name.
            sprintf(name, "./res/mean_dens_");
            sprintf(s_iter, "%d", step);
            strcpy(end, format.c_str());
            strcat(name, s_iter);
            strcat(name, ".");
            strcat(name, end);

            printf("Writing results to file %s...\n", name);

            pngwriter png(n_x, n_y, 0, name);

            Real density;

            // Get the minimum and maximum field data values.
            Real min = 0.0;
            Real max = (Real)n_dir;

            // Write mean density from lattice to file.
            //
            // Loop over all cells.
            for (int cell = 0; cell < n_cells; ++cell)
            {
                // Get row and column index of the current cell.
                int pos_x = cell % n_x;
                int pos_y = cell / n_x;

                if (mean_density_cpu[cell] > -1.0e-06) {

                    density = mean_density_cpu[cell];
                }
#ifdef DEBUG
                else {

                    printf("ERROR in write_results(): "
                           "Negative mean density detected.");
                    abort();
                }
#endif
                // Calculate rgb code.
                std::vector<Real> rgb_code = rgb(min, max, density);

                // Write color data to file.
                png.plot(pos_x+1, pos_y+1, rgb_code[0], rgb_code[1], rgb_code[2]);
            }

            // Close the result file.
            png.close();
        }

	    if (write_mean_velocity == true) {

            // Create file name.
            sprintf(name, "./res/mean_vel_x_");
            sprintf(s_iter, "%d", step);
            strcpy(end, format.c_str());
            strcat(name, s_iter);
            strcat(name, ".");
            strcat(name, end);

            printf("Writing results to file %s...\n", name);

            pngwriter png(n_x, n_y, 0, name);

            Real x_vel;

            // Get the minimum and maximum field data values.
            Real min = - 1.0 * u;
            Real max =   3.0 * u;

            // Write mean velocity in x direction from lattice to file.
            //
            // Loop over all cells.
            for (int cell = 0; cell < n_cells; ++cell)
            {
                // Get row and column index of the current cell.
                int pos_x = cell % n_x;
                int pos_y = cell / n_x;

                if (mean_density_cpu[cell] > 1.0e-06) {

                    x_vel = mean_momentum_cpu[cell] / mean_density_cpu[cell];

                } else if (fabs(mean_density_cpu[cell]) < 1.0e-06) {

                    x_vel = 0.0;
                }

    #ifdef DEBUG

                else {

                    printf("ERROR in write_results(): "
                           "Negative mean density detected.");
                    abort();
                }

    #endif

                // Calculate rgb code.
                std::vector<Real> rgb_code = rgb(min, max, x_vel);

                // Write color data to file.
                png.plot(pos_x+1, pos_y+1, rgb_code[0], rgb_code[1], rgb_code[2]);
            }

            // Close the result file.
            png.close();
	    }

	// Invalid file extension.
	} else {

        printf("ERROR in write_results(): "
               "Invalid file extension %s.", format.c_str());
        abort();
	}

	printf("...done.\n");
}

// Initializes the lattice gas automaton with single particles at defined nodes.
void Lattice::init_single(const std::vector<int> occupied_nodes) {

    // Initialize the lattice with zeros.
    init_zero();

    // Loop over the nodes to place a particle at.
    for (int n = 0; n < occupied_nodes.size(); ++n) {

        node_state_cpu[occupied_nodes[n]] = 1;
    }
}

// Applies boundary conditions for reflecting boundaries at all edges
// of the rectangular domain.
void Lattice::apply_bc_reflecting(const string bounce_type) {

    // Set the cell type to fluid cells.
    int cell_type = 0;

    // Set all cell types to fluid cells.
    apply_cell_type_all(cell_type);

    // Get the type of the boundary cells according to the bounce type.
    if (bounce_type == "back") {

        cell_type = 1;

    } else if (bounce_type == "forward") {

        cell_type = 2;

    } else {

        printf("ERROR in apply_bc_reflecting(): Invalid bounce type %s.\n", bounce_type.c_str());
        abort();
    }

    // Apply the specified cell type to cells located on the boundaries
    // of the rectangular domain.
    apply_boundary_cell_type_east(cell_type);
    apply_boundary_cell_type_north(cell_type);
    apply_boundary_cell_type_west(cell_type);
    apply_boundary_cell_type_south(cell_type);
}

// Applies the specified cell type to all cells in the domain.
void Lattice::apply_cell_type_all(const int cell_type) {

    // Loop over all cells.
#pragma omp parallel for
    for (unsigned int cell = 0; cell < n_cells; ++cell) {

        // Set the cell type to the specified type.
        cell_type_cpu[cell] = cell_type;
    }
}

// Applies the specified cell type to cells located on the eastern boundary
// of the rectangular domain.
void Lattice::apply_boundary_cell_type_east(const int cell_type) {

    // Loop over the cells located at the eastern boundary of
    // the rectangular domain.
    for (int cell = n_x - 1; cell < n_cells; cell += n_x) {

        // Set the cell type to the specified type.
        cell_type_cpu[cell] = cell_type;
    }
}

// Applies the specified cell type to cells located on the northern boundary
// of the rectangular domain.
void Lattice::apply_boundary_cell_type_north(const int cell_type) {

    // Loop over the cells located at the northern boundary of
    // the rectangular domain.
    for (int cell = n_cells - n_x; cell < n_cells; ++cell) {

        // Set the cell type to the specified type.
        cell_type_cpu[cell] = cell_type;
    }
}

// Applies the specified cell type to cells located on the western boundary
// of the rectangular domain.
void Lattice::apply_boundary_cell_type_west(const int cell_type) {

    // Loop over the cells located at the western boundary of
    // the rectangular domain.
    for (int cell = 0; cell < n_cells; cell += n_x) {

        // Set the cell type to the specified type.
        cell_type_cpu[cell] = cell_type;
    }
}

// Applies the specified cell type to cells located on the southern boundary
// of the rectangular domain.
void Lattice::apply_boundary_cell_type_south(const int cell_type) {

    // Loop over the cells located at the southern boundary of
    // the rectangular domain.
    for (int cell = 0; cell < n_x; ++cell) {

        // Set the cell type to the specified type.
        cell_type_cpu[cell] = cell_type;
    }
}

// Prints information about the lattice object to screen.
void Lattice::print_info() {

    printf("Parameter for test case \"%s\":\n", test_case.c_str());
    printf("\n");
    printf("Reynolds number             Re   = %10.2f\n", Re);
    printf("Mach number                 Ma   = %10.2f\n", Ma_s);
    printf("Density                     rho  = %10.2f\n", rho);
    printf("Velocity                    u    = %10.2f\n", u);
    printf("Sound speed                 c    = %10.2f\n", c);
    printf("Scaled sound speed          c_s  = %10.2f\n", c_s);
    printf("Viscosity                   nu   = %10.2f\n", nu);
    printf("Scaled viscosity            nu_s = %10.2f\n", nu_s);
    printf("Galilean breaking factor    g    = %10.2f\n", g);
    printf("\n");
    printf("Number of cells in x direction: %d\n", n_x);
    printf("Number of cells in y direction: %d\n", n_y);
    printf("\n");
}

// Copies all data arrays from the host (CPU) to the device (GPU).
void Lattice::copy_data_to_device()
{
	// Only valid for CUDA parallelized lattice gas automatons.
	// Implemented in CUDA_Lattice.
}

// Copies all data arrays from the device (GPU) back to the host (CPU).
void Lattice::copy_data_from_device()
{
	// Only valid for CUDA parallelized lattice gas automatons.
	// Implemented in CUDA_Lattice.
}

// Computes the number of particles to revert in the context of body force
// in order to accelerate the flow.
int Lattice::get_initial_forcing()
{
	// int equilibrium_forcing = get_equilibrium_forcing();

	return (int)(0.01 * n_cells);
}

// Computes the number of particles to revert in the context of body force.
// in order to compensate boundary layer shear force.
int Lattice::get_equilibrium_forcing()
{
    Real forcing = (8.0 * nu_s * Ma_s * c_s) / pow((Real)n_y, 2.0);

	return ceil(0.5 * n_cells * forcing);
}

// Initializes the lattice gas automaton with some random distributed particles
// in the center area of the domain.
void Lattice::init_diffusion()
{
    // Define the position and size of the center area.
    int  center_x = n_x / 2;
    int  center_y = n_y / 2;
    Real diameter = n_y / 4;

	// Loop over all cells.
#pragma omp parallel for
	for (int cell = 0; cell < n_cells; ++cell) {

		// Get the x and y position of the current cell.
		int pos_x = cell % n_x;
		int pos_y = cell / n_x;

        Real dist = sqrt(pow((pos_x - center_x), 2.0) + pow((pos_y - center_y), 2.0));

	    // Check weather the cell is a fluid cell in the center area of the domain.
	    if (cell_type_cpu[cell] == 0 &&
	    	dist                < (diameter / 2.0))
	    {
            // Loop over all nodes in the fluid cell.
            for (int dir = 0; dir < n_dir; ++dir) {

                // Set random states for the nodes in the fluid cell.
            	node_state_cpu[cell + dir * n_cells] =
            	        (random_uniform() > (1.0 - (1.0 / n_dir)));
            }
	    }
	}
}

// Initializes the lattice gas automaton with some random distributed particles
// in one quarter of a rectangular tank.
void Lattice::init_sloshing()
{
	// Define the width of the water column.
	const int width = n_x / 4;

	// Loop over all cells.
#pragma omp parallel for
	for (int cell = 0; cell < n_cells; ++cell) {

		// Get the x position of the current cell.
		int pos_x = cell % n_x;

	    // Check weather the cell is a fluid cell in the center area of the domain.
	    if (cell_type_cpu[cell] == 0 &&
	    	pos_x                < width )
	    {
            // Loop over all nodes in the fluid cell.
            for (int dir = 0; dir < n_dir; ++dir) {

                // Set random states for the nodes in the fluid cell.
            	node_state_cpu[cell + dir * n_cells] =
            	        (random_uniform() > (1.0 - (1.0 / n_dir)));
            }
	    }
	}
}
