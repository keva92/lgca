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

#include "karman_viewer.h"
#include "ui_karman_viewer.h"

#include "utils.h"
#include "lattice.h"
#include "omp_lattice.h"
#include "cu_lattice.h"
#include "lgca_io_vti.h"
#include "lgca_jet.h"

#include <tbb/task_group.h>

#include <QVTKOpenGLWidget.h>

namespace lgca {

KarmanView::KarmanView(QWidget *parent) :
    QMainWindow(parent),
    m_ui(new Ui::KarmanView),
    m_steps(0)
{
    m_ui->setupUi(this);

    // Print startup message
    print_startup_message();

    // Create a lattice gas cellular automaton object
    m_lattice = new OMP_Lattice<MODEL>(/*case=*/"karman", m_Re, m_Ma, CG_RADIUS);

    // Apply boundary conditions
    m_lattice->apply_bc_karman_vortex_street();

    // Initialize the lattice gas automaton with particles
    m_lattice->init_random();
    m_num_particles = m_lattice->get_n_particles();

    // Accelerate the flow
    m_forcing = m_lattice->get_initial_forcing();
    m_lattice->copy_data_to_output_buffer();
    m_lattice->post_process();
    m_mean_velocity = m_lattice->get_mean_velocity();
    while (m_mean_velocity[0] < m_lattice->u()) {

        // Apply a body force to the particles
        m_lattice->apply_body_force(m_forcing);
        m_lattice->copy_data_to_output_buffer();
        m_lattice->post_process();
        m_mean_velocity = m_lattice->get_mean_velocity();
    }

    // Calculate the number of particles to revert in the context of body force in order to
    // keep the flow up
    m_forcing = m_lattice->get_equilibrium_forcing();

    // Necessary to set up on-line visualization
    m_lattice->copy_data_to_output_buffer();
    m_lattice->post_process();

    // Setup visualization pipeline
    this->setup_visual();

    // Setup UI
    this->setup_ui();

    // Connect widgets
    connect(m_ui->playButton,              SIGNAL(clicked()), this, SLOT(run()));
    connect(m_ui->pauseButton,             SIGNAL(clicked()), this, SLOT(stop()));
    connect(m_ui->rescaleButton,           SIGNAL(clicked()), this, SLOT(rescale()));
    connect(m_ui->cellDensityRadioButton,  SIGNAL(clicked()), this, SLOT(view_cell_density()));
    connect(m_ui->cellMomentumRadioButton, SIGNAL(clicked()), this, SLOT(view_cell_velocity()));
    connect(m_ui->meanDensityRadioButton,  SIGNAL(clicked()), this, SLOT(view_mean_density()));
    connect(m_ui->meanMomentumRadioButton, SIGNAL(clicked()), this, SLOT(view_mean_velocity()));
}

KarmanView::~KarmanView()
{
    m_geom_filter   ->Delete();
    m_mapper        ->Delete();
    m_actor         ->Delete();
    m_ren           ->Delete();
    m_scalar_bar    ->Delete();
    m_scalar_bar_txt->Delete();
    m_lut           ->Delete();
    m_ren_win       ->Delete();
    m_png_filter    ->Delete();
    m_png_writer    ->Delete();

    delete m_vti_io_handler;
    delete m_lattice;
    delete m_ui;
}

void KarmanView::run()
{
    tbb::task_group task_group;

    // Simulation
    task_group.run([&]{

        // Get current mean velocity in x and y direction
        m_mean_velocity = m_lattice->get_mean_velocity();

        // Print current mean velocity in x and y direction
        m_ui->velLineEdit->setText(QStringLiteral("[%1, %2]").arg(
                                       m_mean_velocity[0], /*width=*/5, 'f', /*prec=*/2).arg(
                                       m_mean_velocity[1], /*width=*/5, 'f', /*prec=*/2));
        m_ui->reLineEdit->setText(QString::number(m_lattice->dim_y() / 3 * m_mean_velocity[0] / m_lattice->nu_s(), 'f', /*prec=*/2));
        m_ui->maLineEdit->setText(QString::number(m_mean_velocity[0] / m_lattice->c_s()                          , 'f', /*prec=*/2));

        auto sim_start = steady_clock::now();

#pragma unroll
        for (int s = 0; s < PP_INTERVAL; ++s) {

            // Perform the collision and propagation step on the lattice gas automaton
            m_lattice->collide_and_propagate();
            m_lattice->apply_body_force(m_forcing);
            m_steps++;
        }

        // Print current simulation performance
        auto sim_end = steady_clock::now();
        auto sim_time = std::chrono::duration_cast<duration<double>>(sim_end - sim_start).count();
        m_mnups = (int)((m_lattice->num_cells() * PP_INTERVAL) / (sim_time * 1.0e06));
        m_ui->mnupsLineEdit  ->setText(QString::number(m_mnups));
        m_ui->simTimeLineEdit->setText(QString::number(sim_time, 'f', /*prec=*/2));
        m_ui->stepsLineEdit  ->setText(QString::number(m_steps));

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
        m_ui->ppTimeLineEdit->setText(QString::number(pp_time, 'f', /*prec=*/2));

        // Update image data object
        m_vti_io_handler->update();

        // Update render window
        m_ui->qvtkWidget->GetRenderWindow()->Render();

        // Write image data to file
        if (m_ui->recordButton->isChecked()) {

            if (OUTPUT_FORMAT == "vti") {

                m_vti_io_handler->write(m_steps, OUTPUT_DIR);

            } else if (OUTPUT_FORMAT == "png") {

                std::ostringstream filename;
                filename << OUTPUT_DIR << "res_" << m_steps << ".png";
                m_png_writer->SetFileName(filename.str().c_str());
                m_png_filter->Modified();
                m_png_writer->Write();
            }
        }
    });

    if (!m_ui->pauseButton->isChecked()) QTimer::singleShot(0, this, SLOT(run()));
}

void KarmanView::stop()
{
    // Get the number of particles in the lattice
    unsigned int num_particles_end = m_lattice->get_n_particles();

    // Check weather the number of particles has changed
    if ((num_particles_end - m_num_particles) == 0) {

        fprintf(stderr, "Error check PASSED: There is no difference in the number of particles.\n");

    } else if ((num_particles_end - m_num_particles) != 0) {

        fprintf(stderr, "Error check FAILED: There is a difference in the number of particles of %zu.\n",
                num_particles_end - m_num_particles);
    }
}

void KarmanView::rescale()
{
    double scalarRange[2];
    m_geom_filter->GetOutput()->GetScalarRange(scalarRange);
    m_mapper->SetScalarRange(scalarRange);
}

void KarmanView::view_cell_density()
{
    m_vti_io_handler->set_scalars("Cell density");

    m_geom_filter->SetInputData(m_vti_io_handler->cell_image());
    m_geom_filter->Update();
    m_mapper->SetArrayName("Cell density");
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());

    this->rescale();
    m_ren->ResetCamera();
}

void KarmanView::view_cell_velocity()
{
    m_vti_io_handler->set_scalars("Cell momentum");

    m_geom_filter->SetInputData(m_vti_io_handler->cell_image());
    m_geom_filter->Update();
    m_mapper->SetArrayName("Cell momentum");
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());

    this->rescale();
    m_ren->ResetCamera();
}

void KarmanView::view_mean_density()
{
    m_vti_io_handler->set_scalars("Mean density");

    m_geom_filter->SetInputData(m_vti_io_handler->mean_image());
    m_geom_filter->Update();
    m_mapper->SetArrayName("Mean density");
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());

    this->rescale();
    m_ren->ResetCamera();
}

void KarmanView::view_mean_velocity()
{
    m_vti_io_handler->set_scalars("Mean momentum");

    m_geom_filter->SetInputData(m_vti_io_handler->mean_image());
    m_geom_filter->Update();
    m_mapper->SetArrayName("Mean momentum");
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());

    this->rescale();
    m_ren->ResetCamera();
}

void KarmanView::setup_visual()
{
    m_vti_io_handler = new IoVti<MODEL>(m_lattice, "Mean momentum");

    // Instantiations
    m_geom_filter    = vtkImageDataGeometryFilter::New();
    m_mapper         = vtkPolyDataMapper::New();
    m_actor          = vtkActor::New();
    m_ren            = vtkRenderer::New();
    m_scalar_bar     = vtkScalarBarActor::New();
    m_scalar_bar_txt = vtkTextProperty::New();
    m_lut            = vtkLookupTable::New();
    m_ren_win        = vtkRenderWindow::New();
    m_png_filter     = vtkWindowToImageFilter::New();
    m_png_writer     = vtkPNGWriter::New();

    // Basic pipeline
    m_geom_filter->SetInputData(m_vti_io_handler->mean_image());
    m_geom_filter->Update();
    m_mapper->SetInputConnection(m_geom_filter->GetOutputPort());
    m_mapper->SetArrayName("Mean momentum");
    double scalarRange[2];
    m_geom_filter->GetOutput()->GetScalarRange(scalarRange);
    m_mapper->SetScalarRange(scalarRange);
    m_actor->SetMapper(m_mapper);
    m_ren->AddActor(m_actor);
    m_ren->SetBackground (0.0, 0.0, 0.0);
    m_ren->SetBackground2(0.0, 0.0, 0.0);
    m_ren->GradientBackgroundOn();
    m_ren_win->AddRenderer(m_ren);
    m_ren_win->SetAlphaBitPlanes(1); // Enable usage of alpha channel

    // Scalar bar
    m_scalar_bar_txt->SetFontFamilyToArial();
    m_scalar_bar_txt->SetFontSize(16);
    m_scalar_bar_txt->BoldOff();
    m_scalar_bar_txt->ItalicOff();
    m_scalar_bar_txt->ShadowOff();
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());
    m_scalar_bar->SetNumberOfLabels(5);
    m_scalar_bar->SetTitleTextProperty     (m_scalar_bar_txt);
    m_scalar_bar->SetLabelTextProperty     (m_scalar_bar_txt);
    m_scalar_bar->SetAnnotationTextProperty(m_scalar_bar_txt);
    m_scalar_bar->UnconstrainedFontSizeOn();
    m_scalar_bar->SetBarRatio(0.2);
    m_ren->AddActor2D(m_scalar_bar);

    // Lookup table
    m_lut->SetNumberOfColors(256);
    for (int i = 0; i < m_lut->GetNumberOfColors(); ++i) {
        m_lut->SetTableValue(i, COLORMAP_JET[i][0], COLORMAP_JET[i][1], COLORMAP_JET[i][2], /*a=*/1.0);
    }
    m_lut->SetTableRange(m_mapper->GetScalarRange());
//    m_lut->SetHueRange       (2.0/3.0, 0.0); // Blue to red rainbow
//    m_lut->SetSaturationRange(1.0, 1.0);
//    m_lut->SetValueRange     (1.0, 1.0);
    m_lut->SetVectorModeToMagnitude(); // Does this have any effect?
    m_lut->Build();
    m_mapper    ->SetLookupTable(m_lut);
    m_scalar_bar->SetLookupTable(m_lut);

    // Screenshot
    m_png_filter->SetInput(m_ren_win);
    m_png_filter->SetInputBufferTypeToRGBA();   // Also record the alpha (transparency) channel
    m_png_filter->ReadFrontBufferOff();         // Read from the back buffer

    m_ren->ResetCamera();
}

void KarmanView::setup_ui()
{
    m_ui->qvtkWidget->SetRenderWindow(m_ren_win);
    m_png_filter->Update();
    m_png_writer->SetInputConnection(m_png_filter->GetOutputPort());
    m_ui->qvtkWidget->show();

    m_ui  ->playButton->setIcon(m_ui  ->playButton->style()->standardIcon(QStyle::SP_MediaPlay));
    m_ui ->pauseButton->setIcon(m_ui ->pauseButton->style()->standardIcon(QStyle::SP_MediaPause));
    m_ui->recordButton->setIcon(m_ui->recordButton->style()->standardIcon(QStyle::SP_DialogSaveButton));
    m_ui ->pauseButton->setCheckable(true);
    m_ui->recordButton->setCheckable(true);

    // Set defaults
    m_mean_velocity = m_lattice->get_mean_velocity();
    m_ui->stepsLineEdit       ->setText(QString::number(0));
    m_ui->mnupsLineEdit       ->setText(QString::number(0));
    m_ui->simTimeLineEdit     ->setText(QString::number(0.0, 'f', /*prec=*/2));
    m_ui->ppTimeLineEdit      ->setText(QString::number(0.0, 'f', /*prec=*/2));
    m_ui->velLineEdit         ->setText(QStringLiteral("[%1, %2]").arg(
                                   m_mean_velocity[0], /*width=*/5, 'f', /*prec=*/2).arg(
                                   m_mean_velocity[1], /*width=*/5, 'f', /*prec=*/2));
    m_ui->reLineEdit          ->setText(QString::number(m_lattice->dim_y() / 3 * m_mean_velocity[0] / m_lattice->nu_s(), 'f', /*prec=*/2));
    m_ui->maLineEdit          ->setText(QString::number(m_mean_velocity[0] / m_lattice->c_s()                          , 'f', /*prec=*/2));
    m_ui->numCellsLineEdit    ->setText(QStringLiteral("%1 x %2").arg(m_lattice->dim_x()).arg(m_lattice->dim_y()));
    m_ui->numParticlesLineEdit->setText(QString::number(m_lattice->get_n_particles()));

    m_ui->meanMomentumRadioButton->setChecked(true);
}

} // namespace lgca
