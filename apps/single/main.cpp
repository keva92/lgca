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

#include "cu_lattice.h"
#include "omp_lattice.h"
#include "lgca_io_vti.h"

#include <vtkNew.h>
#include <vtkImageData.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkImageMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkLookupTable.h>
#include <vtkRenderWindow.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkImageViewer.h>

#include <tbb/task_group.h>

#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLWidget.h>

// Main function (CPU/host code)
int main(int argc, char **argv) {

    // Define some variables.
    const int    dim = 2;                // Dimension of the problem (only 2D supported so far)
          Real   Re;                     // Reynolds number
          Real   Ma;                     // Mach number
          int    n_dir;                  // Number of lattice directions
          int    s_max;                  // Number of simulated time steps
          int    coarse_graining_radius; // Coarse graining radius
          int    write_steps;            // Number of steps after which the post-processed results are visualized or written to file
          int    body_force_steps;       // Number of steps after which body force is applied to the particles
          int    body_force_intensity;   // Intensity of the body force
          int    device;                 // ID of the CUDA device to use
          int    max_block_size;         // Maximum CUDA block size in x direction
          string parallel_type;          // Parallelization type (CUDA, OMP)
          string output_format;          // Output format (live or vti)

    // Get values from the command line
    get_vals_from_cmd(argc, argv,
                      &Re, &Ma,
                      &n_dir,
                      &s_max,
                      &coarse_graining_radius,
                      &write_steps,
                      &body_force_steps, &body_force_intensity,
                      &device,
                      &max_block_size,
                      &parallel_type,
                      &output_format);

    srand48(time(NULL));

    // Create some time measurement instances
    Timer *time_measure = new Timer();

    // Print startup message
    print_startup_message();

    // Create a lattice gas cellular automaton object (including allocation of memory on host and device)
    printf("Create lattice gas automaton object and allocate memory...\n");
    Lattice* lattice;
    if (parallel_type == "CUDA") {

//        lattice = new CUDA_Lattice(/*case=*/"collision", Re, Ma, n_dir, coarse_graining_radius, device);

    } else if (parallel_type == "OMP") {

        lattice = new OMP_Lattice(/*case=*/"collision", Re, Ma, n_dir, coarse_graining_radius);

    } else {

        printf("ERROR in main(): Invalid parallelization type %s.\n", parallel_type.c_str());
        abort();
    }
    printf("...done.\n\n");

    // Apply boundary conditions -------------------------------------------------------------------
    printf("Apply boundary conditions...\n");

    lattice->apply_bc_pipe();

    printf("...done.\n\n");


    // Initialize the lattice gas automaton with particles -----------------------------------------
    printf("Initialize the lattice gas automaton...\n");

    lattice->init_single_collision();

    // Necessary to set up on-line visualization
    lattice->copy_data_to_output_buffer();
    lattice->post_process();

    printf("...done.\n\n");

    // Get the number of particles in the lattice
    unsigned int n_particles_start = lattice->get_n_particles();

    if (parallel_type == "CUDA") {

        // Go through all available devices and print their properties to the screen
		device_query();

		// Copy lattice to GPU.
		printf("Copy lattice to device...\n");
		lattice->copy_data_to_device();
		printf("...done.\n\n");
    }

    // Set (proper) parallelization parameters
	lattice->setup_parallel();

    // Initialize timer.
    init_timer(time_measure);


    // Setup on-line visualization -----------------------------------------------------------------

    // Needed to ensure appropriate OpenGL context is created for VTK rendering
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());

    // QT Stuff
    QApplication app(argc, argv);
    QApplication::setStyle("fusion");

    lgca::IoVti vtiIoHandler(lattice, "Cell density");

    vtkNew<vtkImageDataGeometryFilter> geomFilter;
    geomFilter->SetInputData(vtiIoHandler.image());
    geomFilter->Update();
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(geomFilter->GetOutputPort());
    mapper->SetArrayName("Cell density");
    double scalarRange[2];
    geomFilter->GetOutput()->GetScalarRange(scalarRange);
    mapper->SetScalarRange(scalarRange);
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper.GetPointer());
    vtkNew<vtkRenderer> ren;
    ren->AddActor(actor.GetPointer());
    ren->SetBackground(1.0, 1.0, 1.0);
    ren->SetBackground2(0.2, 0.3, 0.5);
    ren->GradientBackgroundOn();
    vtkNew<vtkScalarBarActor> scalarBar;
    vtkNew<vtkTextProperty> text;
    text->SetFontFamilyToArial();
    text->SetFontSize(10);
    text->BoldOff();
    text->ItalicOff();
    text->ShadowOff();
    scalarBar->SetTitle(mapper->GetArrayName());
    scalarBar->SetNumberOfLabels(4);
    scalarBar->SetTitleTextProperty(text.GetPointer());
    scalarBar->SetLabelTextProperty(text.GetPointer());
    scalarBar->SetBarRatio(0.2);
    ren->AddActor2D(scalarBar.GetPointer());
    vtkNew<vtkLookupTable> lut;
    lut->SetTableRange(mapper->GetScalarRange());
    lut->SetHueRange(2.0/3.0, 0.0); // Blue to red rainbow
    lut->SetSaturationRange(1.0, 1.0);
    lut->SetValueRange(1.0, 1.0);
    lut->Build();
    mapper->SetLookupTable(lut.GetPointer());
    scalarBar->SetLookupTable(lut.GetPointer());

    vtkNew<vtkRenderWindow> renWin;
    renWin->AddRenderer(ren.GetPointer());
    renWin->SetSize(1024, 576); // WXGA
    ren->ResetCamera();

    // TODO Write GUI
//    QVTKOpenGLWidget widget;
//    widget.resize(1024, 576); // WXGA
//    vtkNew<vtkGenericOpenGLRenderWindow> glRenWin;
//    glRenWin->AddRenderer(ren.GetPointer());
//    widget.SetRenderWindow(glRenWin.GetPointer());
//    widget.show();
//    app.exec();


    // Time loop -----------------------------------------------------------------------------------

    tbb::task_group task_group;

    printf("Starting calculation...\n");

    auto total_start = steady_clock::now();

    // Loop over bunches of simulation time steps
    for (int step = 0; step <= s_max/write_steps; ++step) {

        // Simulation
        task_group.run([&]{

            auto sim_start = steady_clock::now();

            for (int s = 0; s < write_steps; ++s) {

                // Perform the collision and propagation step on the lattice gas automaton
                lattice->collide_and_propagate(step);
            }

            // Print current simulation performance
            auto sim_end = steady_clock::now();
            auto sim_time = std::chrono::duration_cast<duration<double>>(sim_end - sim_start).count();
            fprintf(stderr, "Current MNUPS: %d\n", (int)((lattice->get_n_x() * lattice->get_n_y() * write_steps)
                                                      / (sim_time * 1.0e06)));

            // Copy results back to host
            if (parallel_type == "CUDA") lattice->copy_data_from_device();

            // Copy results to a temporary buffer for post-processing and visualization
            lattice->copy_data_to_output_buffer();
        });

        // Visualization
        task_group.run_and_wait([&]{

            // Compute quantities of interest as a post-processing procedure
            lattice->post_process();

            // Update image data object
            vtiIoHandler.update();

            // Update render window
            renWin->Render();
            sleep(1);
        });

    } // for steps

    auto total_end = steady_clock::now();

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

    auto total_time = std::chrono::duration_cast<duration<double>>(total_end - total_start).count();

    printf("Total elapsed time: %e s for %d simulation steps.\n", total_time, s_max);
    printf("Average MNUPS: %d\n", (int)((lattice->get_n_x() * lattice->get_n_y() * s_max)
                                     / (total_time * 1.0e06)));
    printf("\n");

    printf("Terminate program...\n");
    delete lattice;
    printf("...done.\n");

    // Exit program.
    return EXIT_SUCCESS;
}
