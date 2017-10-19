/*
 * main.cpp
 *
 *  Created on: Nov 26, 2015
 *      Author: Kerstin Vater
 * Description: This is a lattice gas cellular automaton in two dimensions.
 *
 */

// C++ include files
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <malloc.h>
#include <math.h>
#include <vector>

// User-defined include files
#include "cu_lattice.h"
#include "omp_lattice.h"

// Namespaces
using namespace std;

// Define real number precision to use.
typedef float real;

// Main function (CPU/host code)
int main(int argc, char **argv) {

    // Define some variables.
    const int    dim = 2;                // Dimension of the problem.
          string test_case;              // Test case
                                         // (pipe flow,
          	  	  	  	  	  	  	  	 // box,
                                         // Karman vortex street,
          	  	  	  	  	  	  	  	 // single collision,
          	  	  	  	  	  	  	  	 // diffusion,
          	  	  	  	  	  	  	  	 // sloshing,
          	  	  	  	  	  	  	  	 // hourglass).
          real   Re;                     // Reynolds number.
          real   Ma;                     // Mach number.
          int    n_dir;                  // Number of lattice directions.
          int    s_max;                  // Number of simulated time steps.
          int    coarse_graining_radius; // Coarse graining radius.
          int    write_steps;            // Number of steps after which the post-processed results are written to a file.
          int    body_force_steps;       // Number of steps after which a body force is applied to the particles.
          int    body_force_intensity;   // Intensity of the body force.
          int    device;                 // Number of the device to use.
          int    max_block_size;         // Maximum block size in x direction.
          string parallel_type;          // Parallelization type (CUDA, openMP).
          string output_format;          // Output format (png or vti).

    // Get values from the command line.
    get_vals_from_cmd(argc, argv,
                      &test_case,
                      &Re, &Ma,
                      &n_dir,
                      &s_max,
                      &coarse_graining_radius,
                      &write_steps,
                      &body_force_steps, &body_force_intensity,
                      &device,
                      &max_block_size,
                      &parallel_type,
                      & output_format);

    srand48(time(NULL));

    // Create a new time measurement instance.
    Timer *time_measure = new Timer();

    // Print startup message.
    print_startup_message();

    // Create a lattice gas cellular automaton object
    // (including allocation of memory on host and device).
    printf("Create lattice gas automaton object and allocate memory...\n");
    Lattice* lattice;
    if (parallel_type == "CUDA") {

    	lattice = new CUDA_Lattice(test_case, Re, Ma, n_dir, coarse_graining_radius, device);

    } else if (parallel_type == "openMP") {

    	lattice = new OMP_Lattice(test_case, Re, Ma, n_dir, coarse_graining_radius);

    } else {

        printf("ERROR in main(): Invalid parallelization type %s.\n", parallel_type.c_str());
        abort();
    }
    printf("...done.\n\n");

    // Apply boundary conditions.
    printf("Apply boundary conditions...\n");

    if (test_case == "pipe") {

        lattice->apply_bc_pipe();

    } else if (test_case == "karman") {

        lattice->apply_bc_karman_vortex_street();

    } else if (test_case == "box" || test_case == "diffusion" || test_case == "sloshing") {

    	lattice->apply_bc_reflecting("back");

    } else if (test_case == "collision") {

        lattice->apply_bc_pipe();

    } else if (test_case == "periodic") {

    	lattice->apply_bc_periodic();

    } else if (test_case == "hourglass") {

    	// lattice->apply_bc_hourglass();

    } else {

        printf("ERROR in main(): Invalid test case %s.\n", test_case.c_str());
        abort();
    }

    printf("...done.\n\n");

    // Initialize the lattice gas automaton with particles.
    printf("Initialize the lattice gas automaton...\n");

    if ((test_case == "pipe") || (test_case == "karman") || (test_case == "box") || (test_case == "periodic")) {

        lattice->init_random();

    } else if (test_case == "collision") {

        lattice->init_single_collision();

    } else if (test_case == "diffusion") {

    	lattice->init_diffusion();

    } else if (test_case == "sloshing") {

    	lattice->init_sloshing();

    } else if (test_case == "hourglass") {

    	// lattice->init_hourglass();

    } else {

        printf("ERROR in main(): Invalid test case %s.\n", test_case.c_str());
        abort();
    }

    printf("...done.\n\n");

    // Get the number of particles in the lattice.
    unsigned int n_particles_start = lattice->get_n_particles();

    if (parallel_type == "CUDA") {

		// Go through all available devices and print their properties to the screen.
		device_query();

		// Copy lattice to GPU.
		printf("Copy lattice to device...\n");
		lattice->copy_data_to_device();
		printf("...done.\n\n");
    }

	// Set (proper) parallelization parameters.
	lattice->setup_parallel();

    // Initialize timer.
    init_timer(time_measure);

    // Define variables for performance measure.
    real time_start;
    real time_end;

    vector<real> mean_velocity;

    // Calculate the number of particles to revert in the context of body force
    // in order to accelerate the flow.
    int forcing = lattice->get_initial_forcing();

    printf("Starting calculation...\n");

    // Loop over the simulations steps.
    for (unsigned int step = 0; step <= s_max; ++step) {

        if (step % body_force_steps == 0 ||
            step % write_steps      == 0)
        {
            // Compute quantities of interest as a post-processing procedure.
            lattice->post_process();

            // Copy results back to host.
            if (parallel_type == "CUDA") lattice->copy_data_from_device();

            // Get current mean velocity in x and y direction.
            mean_velocity = lattice->get_mean_velocity();
        }

        if (step % write_steps == 0) {

            printf("\n");
            printf("Executing step %d...\n", step);

            // Print current mean velocity in x and y direction.
            printf("Current mean velocity: (%6.4f, %6.4f)\n", mean_velocity[0], mean_velocity[1]);

            // Get current performance.
            if (step > 0) {

                // Get time.
                time_end = get_timer(time_measure);

                // Print current performance.
                printf("Current MNUPS: %d\n", (int) ((lattice->get_n_x() * lattice->get_n_y() * write_steps) / ((time_end - time_start) * 1.0e06)));
            }

            // Write results of current time step to file.
            lattice->write_results(step, output_format);

#ifdef DEBUG

            // Check weather the number of particles in the lattice has changed.
            printf("There are %d particles in the lattice.\n", lattice->get_n_particles());

            // Print the lattice to the screen.
            lattice->print();
#endif

            // Get time.
            time_start = get_timer(time_measure);
        }

        // Call the CUDA kernel which performs the collision and propagation step
        // on the lattice gas automaton.
        lattice->collide_and_propagate(step);

        if ((body_force_intensity > 1.0e-06) &&
        	(mean_velocity[0] < lattice->get_u()) &&
        	(test_case == "pipe" || test_case == "karman")) {

        	// Reduce the forcing once the flow has been accelerated strong enough.
        	if (mean_velocity[0] > 0.9 * lattice->get_u()) forcing = lattice->get_equilibrium_forcing();

            // Apply a body force to the particles.
            lattice->apply_body_force(forcing);
        }

        if ((body_force_intensity > 1.0e-06) &&
        	(test_case == "sloshing" || test_case == "hourglass")) {

            // Apply a body force to the particles.
            lattice->apply_body_force(body_force_intensity);
        }

    }
    printf("\n");

    // Get the number of particles in the lattice.
    unsigned int n_particles_end = lattice->get_n_particles();

    // Check weather the number of particles has changed.
    if ((n_particles_end - n_particles_start) == 0) {

        printf("Error check PASSED: There is no difference in the number of particles.\n");

    } else if ((n_particles_end - n_particles_start) != 0) {

        printf("Error check FAILED: There is a difference in the number of particles of %d.\n", n_particles_end - n_particles_start);
    }
    printf("\n");

    // Stop time measurement.
    real total_elapsed_time = get_timer(time_measure);

    printf("Total elapsed time: %e s for %d simulation steps.\n", total_elapsed_time, s_max);
    printf("Average MNUPS: %d\n", (int) ((lattice->get_n_x() * lattice->get_n_y() * s_max) / (total_elapsed_time * 1.0e06)));
    printf("\n");

    printf("Terminate program...\n");
    delete lattice;
    printf("...done.\n");

    // Exit program.
    return EXIT_SUCCESS;
}
