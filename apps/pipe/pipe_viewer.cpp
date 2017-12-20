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

#include "pipe_viewer.h"
#include "ui_pipe_viewer.h"

#include "utils.h"
#include "cuda_utils.cuh"
#include "lattice.h"
#include "omp_lattice.h"
#include "cu_lattice.h"
#include "lgca_io_vti.h"

#include <tbb/task_group.h>

#include <QVTKOpenGLWidget.h>

namespace lgca {

PipeView::PipeView(QWidget *parent) :
    QMainWindow(parent),
    m_ui(new Ui::PipeView)
{
    m_ui->setupUi(this);

    // Define some variables
    Real   Re                     = 80.0;     // Reynolds number
    Real   Ma                     = 0.2;      // Mach number
    int    coarse_graining_radius = 15;       // Coarse graining radius

    srand48(time(NULL));

    // Print startup message
    print_startup_message();

    // Create a lattice gas cellular automaton object
    m_lattice = new OMP_Lattice<NUM_DIR>(/*case=*/"pipe", Re, Ma, coarse_graining_radius);

    // Apply boundary conditions
    m_lattice->apply_bc_pipe();

    // Initialize the lattice gas automaton with particles
    m_lattice->init_random();
    m_num_particles = m_lattice->get_n_particles();

    // Necessary to set up on-line visualization
    m_lattice->copy_data_to_output_buffer();
    m_lattice->post_process();

    // Calculate the number of particles to revert in the context of body force in order to
    // accelerate the flow.
    m_forcing = m_lattice->get_initial_forcing();

    // Set (proper) parallelization parameters
    m_lattice->setup_parallel();

    m_vti_io_handler = new IoVti<NUM_DIR>(m_lattice, "Mean momentum");

    m_geom_filter   = vtkImageDataGeometryFilter::New();
    m_mapper        = vtkPolyDataMapper::New();
    m_actor         = vtkActor::New();
    m_ren           = vtkRenderer::New();
    m_scalar_bar    = vtkScalarBarActor::New();
    m_scalar_bar_txt = vtkTextProperty::New();
    m_lut           = vtkLookupTable::New();
    m_ren_win       = vtkRenderWindow::New();

    m_geom_filter->SetInputData(m_vti_io_handler->mean_image());
    m_geom_filter->Update();
    m_mapper->SetInputConnection(m_geom_filter->GetOutputPort());
    m_mapper->SetArrayName("Mean momentum");
    double scalarRange[2];
    m_geom_filter->GetOutput()->GetScalarRange(scalarRange);
    m_mapper->SetScalarRange(scalarRange);
    m_actor->SetMapper(m_mapper);
    m_ren->AddActor(m_actor);
    m_ren->SetBackground(1.0, 1.0, 1.0);
    m_ren->SetBackground2(0.2, 0.3, 0.5);
    m_ren->GradientBackgroundOn();
    m_scalar_bar_txt->SetFontFamilyToArial();
    m_scalar_bar_txt->SetFontSize(10);
    m_scalar_bar_txt->BoldOff();
    m_scalar_bar_txt->ItalicOff();
    m_scalar_bar_txt->ShadowOff();
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());
    m_scalar_bar->SetNumberOfLabels(4);
    m_scalar_bar->SetTitleTextProperty(m_scalar_bar_txt);
    m_scalar_bar->SetLabelTextProperty(m_scalar_bar_txt);
    m_scalar_bar->SetBarRatio(0.2);
    m_ren->AddActor2D(m_scalar_bar);
    m_lut->SetTableRange(m_mapper->GetScalarRange());
    m_lut->SetHueRange(2.0/3.0, 0.0); // Blue to red rainbow
    m_lut->SetSaturationRange(1.0, 1.0);
    m_lut->SetValueRange(1.0, 1.0);
    m_lut->Build();
    m_mapper->SetLookupTable(m_lut);
    m_scalar_bar->SetLookupTable(m_lut);
    m_ren->ResetCamera();

    m_ui->qvtkWidget->SetRenderWindow(m_ren_win);
    m_ui->qvtkWidget->GetRenderWindow()->AddRenderer(m_ren);
    m_ui->qvtkWidget->show();

    connect(m_ui->startButton, SIGNAL(clicked()), this, SLOT(run()));
    connect(m_ui->stopButton, SIGNAL(clicked()), this, SLOT(stop()));
    connect(m_ui->rescaleButton, SIGNAL(clicked()), this, SLOT(rescale()));
}

PipeView::~PipeView()
{
    m_geom_filter   ->Delete();
    m_mapper        ->Delete();
    m_actor         ->Delete();
    m_ren           ->Delete();
    m_scalar_bar    ->Delete();
    m_scalar_bar_txt->Delete();
    m_lut           ->Delete();
    m_ren_win       ->Delete();

    delete m_vti_io_handler;
    delete m_lattice;
    delete m_ui;
}

void PipeView::run()
{
    tbb::task_group task_group;

    // Simulation
    task_group.run([&]{

        // Get current mean velocity in x and y direction
        m_mean_velocity = m_lattice->get_mean_velocity();

        if (m_mean_velocity[0] < m_lattice->u()) {

            // Reduce the forcing once the flow has been accelerated strong enough
            if (m_mean_velocity[0] > 0.9 * m_lattice->u())
                m_forcing = m_lattice->get_equilibrium_forcing();

            // Apply a body force to the particles
            m_lattice->apply_body_force(m_forcing);
        }

        auto sim_start = steady_clock::now();

        for (int s = 0; s < m_write_steps; ++s) {

            // Perform the collision and propagation step on the lattice gas automaton
            m_lattice->collide_and_propagate(s);
        }

        // Print current simulation performance
        auto sim_end = steady_clock::now();
        auto sim_time = std::chrono::duration_cast<duration<double>>(sim_end - sim_start).count();
        fprintf(stderr, "Simulation took %f s.\n", sim_time);
        m_mnups = (int)((m_lattice->num_cells() * m_write_steps) / (sim_time * 1.0e06));
        fprintf(stderr, "Current MNUPS: %d\n", m_mnups);

        // Print current mean velocity in x and y direction
        fprintf(stderr, "Current mean velocity: (%6.4f, %6.4f)\n", m_mean_velocity[0], m_mean_velocity[1]);

        // Copy results to a temporary buffer for post-processing and visualization
        m_lattice->copy_data_to_output_buffer();
    });

    // Visualization
    task_group.run_and_wait([&]{


        // Compute quantities of interest as a post-processing procedure
        auto pp_start = steady_clock::now();
        m_lattice->post_process();
        auto pp_end = steady_clock::now();
        auto pp_time = std::chrono::duration_cast<duration<double>>(pp_end - pp_start).count();
        fprintf(stderr, "Postprocessing took %f s.\n", pp_time);

        // Update image data object
        m_vti_io_handler->update();

        // Update render window
        auto ren_start = steady_clock::now();
        m_ui->qvtkWidget->GetRenderWindow()->Render();
        auto ren_end = steady_clock::now();
        auto ren_time = std::chrono::duration_cast<duration<double>>(ren_end - ren_start).count();
        fprintf(stderr, "Rendering took %f s.\n", ren_time);
    });

    QTimer::singleShot(0, this, SLOT(run()));
}

void PipeView::stop()
{
    // Get the number of particles in the lattice.
    unsigned int num_particles_end = m_lattice->get_n_particles();

    // Check weather the number of particles has changed.
    if ((num_particles_end - m_num_particles) == 0) {

        printf("Error check PASSED: There is no difference in the number of particles.\n");

    } else if ((num_particles_end - m_num_particles) != 0) {

        printf("Error check FAILED: There is a difference in the number of particles of %d.\n", num_particles_end - m_num_particles);
    }

    qApp->quit();
}

void PipeView::rescale()
{
    double scalarRange[2];
    m_geom_filter->GetOutput()->GetScalarRange(scalarRange);
    m_mapper->SetScalarRange(scalarRange);
}

} // namespace lgca
