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

#include "lgca_viewer.h"
#include "ui_lgca_viewer.h"

#include "utils.h"
#include "cuda_utils.cuh"
#include "lattice.h"
#include "omp_lattice.h"
#include "cu_lattice.h"
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

#include <QVTKOpenGLWidget.h>

namespace lgca {

LgcaView::LgcaView(QWidget *parent) :
    QMainWindow(parent),
    m_ui(new Ui::LgcaView)
{
    m_ui->setupUi(this);

    // Define some variables
    Real   Re                     = 80.0;     // Reynolds number
    Real   Ma                     = 0.2;      // Mach number
    int    n_dir                  = 6;        // Number of lattice directions
    int    coarse_graining_radius = 15;       // Coarse graining radius

    srand48(time(NULL));

    // Print startup message
    print_startup_message();

    // Create a lattice gas cellular automaton object
    printf("Create lattice gas automaton object and allocate memory...\n");
    m_lattice = new OMP_Lattice(/*case=*/"collision", Re, Ma, n_dir, coarse_graining_radius);
    printf("...done.\n\n");

    // Apply boundary conditions
    printf("Apply boundary conditions...\n");
    m_lattice->apply_bc_pipe();
    printf("...done.\n\n");

    // Initialize the lattice gas automaton with particles
    printf("Initialize the lattice gas automaton...\n");
    m_lattice->init_single_collision();

    // Necessary to set up on-line visualization
    m_lattice->copy_data_to_output_buffer();
    m_lattice->post_process();
    printf("...done.\n\n");

    // Set (proper) parallelization parameters
    m_lattice->setup_parallel();

    m_vtiIoHandler = new IoVti(m_lattice, "Cell density");

    vtkNew<vtkImageDataGeometryFilter> geomFilter;
    geomFilter->SetInputData(m_vtiIoHandler->image());
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
    m_ui->qvtkWidget->SetRenderWindow(renWin.GetPointer());
    m_ui->qvtkWidget->GetRenderWindow()->AddRenderer(ren.GetPointer());
    ren->ResetCamera();
    m_ui->qvtkWidget->show();

    connect(m_ui->startButton, SIGNAL(clicked()), this, SLOT(run()));
}

LgcaView::~LgcaView()
{
    delete m_vtiIoHandler;
    delete m_lattice;
    delete m_ui;
}

void LgcaView::run()
{
    int write_steps = 1; // Number of steps after which the post-processed results are visualized or written to file

    // Get the number of particles in the lattice
    unsigned int n_particles_start = m_lattice->get_n_particles();

    tbb::task_group task_group;

    printf("Starting calculation...\n");

    auto total_start = steady_clock::now();

    // Loop over bunches of simulation time steps
    int step_count = 0;
    for (int step = 0; step < 10; ++step) {

        // Simulation
        task_group.run([&]{

            auto sim_start = steady_clock::now();

            for (int s = 0; s < write_steps; ++s, ++step_count) {

                // Perform the collision and propagation step on the lattice gas automaton
                m_lattice->collide_and_propagate(s);
            }

            // Print current simulation performance
            auto sim_end = steady_clock::now();
            auto sim_time = std::chrono::duration_cast<duration<double>>(sim_end - sim_start).count();
            fprintf(stderr, "Current MNUPS: %d\n", (int)((m_lattice->num_cells() * write_steps)
                                                       / (sim_time * 1.0e06)));

            // Copy results to a temporary buffer for post-processing and visualization
            m_lattice->copy_data_to_output_buffer();
        });

        // Visualization
        task_group.run_and_wait([&]{

            // Compute quantities of interest as a post-processing procedure
            m_lattice->post_process();

            // Update image data object
            m_vtiIoHandler->update();

            // Update render window
            m_ui->qvtkWidget->GetRenderWindow()->Render();
            sleep(1);
        });

    } // for steps

    auto total_end = steady_clock::now();

    printf("\n");

    // Get the number of particles in the lattice.
    unsigned int n_particles_end = m_lattice->get_n_particles();

    // Check weather the number of particles has changed.
    if ((n_particles_end - n_particles_start) == 0) {

        printf("Error check PASSED: There is no difference in the number of particles.\n");

    } else if ((n_particles_end - n_particles_start) != 0) {

        printf("Error check FAILED: There is a difference in the number of particles of %d.\n", n_particles_end - n_particles_start);
    }
    printf("\n");

    auto total_time = std::chrono::duration_cast<duration<double>>(total_end - total_start).count();

    printf("Total elapsed time: %e s for %d simulation steps.\n", total_time, step_count);
    printf("Average MNUPS: %d\n", (int)((m_lattice->num_cells() * step_count)
                                      / (total_time * 1.0e06)));
    printf("\n");

    printf("Terminate program...\n");
    printf("...done.\n");
}

} // namespace lgca
