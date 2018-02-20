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
#include "lattice.h"
#include "omp_lattice.h"
#include "cu_lattice.h"
#include "lgca_io_vti.h"
#include "lgca_jet.h"

#include <tbb/task_group.h>

#include <QVTKOpenGLWidget.h>

namespace lgca {

PipeView::PipeView(QWidget *parent) :
    QMainWindow(parent),
    m_ui(new Ui::PipeView),
    m_steps(0)
{
    m_ui->setupUi(this);

    // Connect widgets
    connect(m_ui->actionSingle,     &QAction::triggered, this, &PipeView::setup_single);
    connect(m_ui->actionDiffusion,  &QAction::triggered, this, &PipeView::setup_diffusion);
    connect(m_ui->actionPipe,       &QAction::triggered, this, &PipeView::setup_pipe);
    connect(m_ui->actionKarman,     &QAction::triggered, this, &PipeView::setup_karman);
    connect(m_ui->actionExit,       &QAction::triggered, this, &PipeView::quit);
    connect(m_ui->actionAbout,      &QAction::triggered, this, &PipeView::show_about);

    connect(m_ui->playButton,              SIGNAL(clicked()),   this, SLOT(run()));
    connect(m_ui->pauseButton,             SIGNAL(clicked()),   this, SLOT(stop()));
    connect(m_ui->rescaleButton,           SIGNAL(clicked()),   this, SLOT(rescale()));
    connect(m_ui->cellDensityRadioButton,  SIGNAL(clicked()),   this, SLOT(view_cell_density()));
    connect(m_ui->cellMomentumRadioButton, SIGNAL(clicked()),   this, SLOT(view_cell_momentum()));
    connect(m_ui->meanDensityRadioButton,  SIGNAL(clicked()),   this, SLOT(view_mean_density()));
    connect(m_ui->meanMomentumRadioButton, SIGNAL(clicked()),   this, SLOT(view_mean_momentum()));
    connect(m_ui->streamlinesRadioButton,  SIGNAL(clicked()),   this, SLOT(view_streamlines()));

    m_ui->centralWidget->hide();
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
    m_png_filter    ->Delete();
    m_png_writer    ->Delete();
    m_stream_filter ->Delete();
    m_streak_filter ->Delete();
    m_path_filter   ->Delete();

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

        // Print current mean velocity in x and y direction
        m_ui->velLineEdit->setText(QStringLiteral("[%1, %2]").arg(
                                       m_mean_velocity[0], /*width=*/5, 'f', /*prec=*/2).arg(
                                       m_mean_velocity[1], /*width=*/5, 'f', /*prec=*/2));
        m_ui->reLineEdit->setText(QString::number(m_lattice->dim_y() * m_mean_velocity[0] / m_lattice->nu_s(), 'f', /*prec=*/2));
        m_ui->maLineEdit->setText(QString::number(m_mean_velocity[0] / m_lattice->c_s()                      , 'f', /*prec=*/2));

        if (m_mean_velocity[0] < m_lattice->u()) {

            // Reduce the forcing once the flow has been accelerated strong enough
            if (m_mean_velocity[0] > 0.9 * m_lattice->u())
                m_forcing = m_lattice->get_equilibrium_forcing();

            // Apply a body force to the particles
            m_lattice->apply_body_force(m_forcing);
        }

        auto sim_start = steady_clock::now();

#pragma unroll
        for (int s = 0; s < PP_INTERVAL; ++s) {

            // Perform the collision and propagation step on the lattice gas automaton
            m_lattice->collide_and_propagate();
            m_steps++;
        }

        // Print current simulation performance
        auto sim_end = steady_clock::now();
        auto sim_time = std::chrono::duration_cast<duration<double>>(sim_end - sim_start).count();
        m_mnups = (int)((m_lattice->num_cells() * PP_INTERVAL) / (sim_time * 1.0e06));
        m_ui->mnupsLineEdit  ->setText(QString::number(m_mnups));
        m_ui->simTimeLineEdit->setText(QString::number(sim_time, 'f', /*prec=*/2));
        m_ui->stepsLineEdit  ->setText(QString::number(m_steps));
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

    // Copy results to a temporary buffer for post-processing and visualization
    m_lattice->copy_data_to_output_buffer();

    if (!m_ui->pauseButton->isChecked()) QTimer::singleShot(0, this, SLOT(run()));
}

void PipeView::stop()
{
    // Get the number of particles in the lattice
    size_t num_particles_end = m_lattice->get_n_particles();

    // Check weather the number of particles has changed
    if ((num_particles_end - m_num_particles) == 0) {

        fprintf(stderr, "Error check PASSED: There is no difference in the number of particles.\n");

    } else if ((num_particles_end - m_num_particles) != 0) {

        fprintf(stderr, "Error check FAILED: There is a difference in the number of particles of %zu.\n",
                num_particles_end - m_num_particles);
    }
}

void PipeView::rescale()
{
    double scalarRange[2];
    m_geom_filter->GetOutput()->GetScalarRange(scalarRange); // TODO
    m_mapper->SetScalarRange(scalarRange);
    m_ui->qvtkWidget->GetRenderWindow()->Render();
}

void PipeView::view_cell_density()
{
    m_vti_io_handler->set_scalars("Cell density");

    m_geom_filter->SetInputData(m_vti_io_handler->cell_image());
    m_geom_filter->Update();
    m_mapper->SetInputConnection(m_geom_filter->GetOutputPort());
    m_mapper->SetArrayName("Cell density");
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());

    this->rescale();
    m_ren->ResetCamera();
    m_ui->qvtkWidget->GetRenderWindow()->Render();
}

void PipeView::view_cell_momentum()
{
    m_vti_io_handler->set_scalars("Cell momentum");

    m_geom_filter->SetInputData(m_vti_io_handler->cell_image());
    m_geom_filter->Update();
    m_mapper->SetInputConnection(m_geom_filter->GetOutputPort());
    m_mapper->SetArrayName("Cell momentum");
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());

    this->rescale();
    m_ren->ResetCamera();
    m_ui->qvtkWidget->GetRenderWindow()->Render();
}

void PipeView::view_mean_density()
{
    m_vti_io_handler->set_scalars("Mean density");

    m_geom_filter->SetInputData(m_vti_io_handler->mean_image());
    m_geom_filter->Update();
    m_mapper->SetInputConnection(m_geom_filter->GetOutputPort());
    m_mapper->SetArrayName("Mean density");
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());

    this->rescale();
    m_ren->ResetCamera();
    m_ui->qvtkWidget->GetRenderWindow()->Render();
}

void PipeView::view_mean_momentum()
{
    m_vti_io_handler->set_scalars("Mean momentum");

    m_geom_filter->SetInputData(m_vti_io_handler->mean_image());
    m_geom_filter->Update();
    m_mapper->SetInputConnection(m_geom_filter->GetOutputPort());
    m_mapper->SetArrayName("Mean momentum");
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());

    this->rescale();
    m_ren->ResetCamera();
    m_ui->qvtkWidget->GetRenderWindow()->Render();
}

void PipeView::view_streamlines()
{
    m_vti_io_handler->set_vectors("Mean momentum");

    m_stream_filter->SetInputData(m_vti_io_handler->mean_image());
    m_stream_filter->Update();
    m_mapper->SetInputConnection(m_stream_filter->GetOutputPort());
    m_mapper->SetArrayName("Mean momentum");
    m_scalar_bar->SetTitle(m_mapper->GetArrayName());

    this->rescale();
    m_ren->ResetCamera();
    m_ui->qvtkWidget->GetRenderWindow()->Render();
}

void PipeView::setup_visual()
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
    m_stream_filter  = vtkStreamTracer::New();
    m_streak_filter  = vtkStreaklineFilter::New();
    m_path_filter    = vtkParticlePathFilter::New();

    // Streamlines
    vtkLineSource* line = vtkLineSource::New();
    line->SetPoint1(0.0, 0.0, 0.0);
    line->SetPoint2(0.0, 10.0, 0.0);
    line->SetResolution(40);
    line->Update();
    m_stream_filter->SetInputData(m_vti_io_handler->mean_image());
    m_stream_filter->SetSourceData(line->GetOutput());
    m_stream_filter->SetMaximumPropagation(100.0);
    m_stream_filter->SetComputeVorticity(true);
    m_stream_filter->Update();
    line->Delete();

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
    m_lut->Build();
    m_mapper    ->SetLookupTable(m_lut);
    m_scalar_bar->SetLookupTable(m_lut);

    // Screenshot
    m_png_filter->SetInput(m_ren_win);
    m_png_filter->SetInputBufferTypeToRGBA();   // Also record the alpha (transparency) channel
    m_png_filter->ReadFrontBufferOff();         // Read from the back buffer

    m_ren->ResetCamera();
}

void PipeView::setup_ui()
{
    m_ui->qvtkWidget->SetRenderWindow(m_ren_win);
    m_png_filter->Update();
    m_png_writer->SetInputConnection(m_png_filter->GetOutputPort());
    m_ui->qvtkWidget->show();

    m_ui  ->playButton->setIcon(m_ui->  playButton->style()->standardIcon(QStyle::SP_MediaPlay));
    m_ui ->pauseButton->setIcon(m_ui-> pauseButton->style()->standardIcon(QStyle::SP_MediaPause));
    m_ui->recordButton->setIcon(m_ui->recordButton->style()->standardIcon(QStyle::SP_DialogSaveButton));
    m_ui->pauseButton->setCheckable(true);
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
    m_ui->reLineEdit          ->setText(QString::number(m_lattice->dim_y() * m_mean_velocity[0] / m_lattice->nu_s(), 'f', /*prec=*/2));
    m_ui->maLineEdit          ->setText(QString::number(m_mean_velocity[0] / m_lattice->c_s()                      , 'f', /*prec=*/2));
    m_ui->numCellsLineEdit    ->setText(QStringLiteral("%1 x %2").arg(m_lattice->dim_x()).arg(m_lattice->dim_y()));
    m_ui->numParticlesLineEdit->setText(QString::number(m_lattice->get_n_particles()));

    m_ui->meanMomentumRadioButton->setChecked(true);
}

void PipeView::setup_single()
{
    // TODO
    cout << "Setup single particle collision case..." << endl;
}

void PipeView::setup_diffusion()
{
    // TODO
    cout << "Setup molecular diffusion case..." << endl;
}

void PipeView::setup_pipe()
{
    // Create a lattice gas cellular automaton object
    m_lattice = new OMP_Lattice<MODEL>(/*case=*/"pipe", m_Re, m_Ma, CG_RADIUS);

    // Apply boundary conditions
    m_lattice->apply_bc_pipe();

    // Initialize the lattice gas automaton with particles
    m_lattice->init_pipe();
    m_num_particles = m_lattice->get_n_particles();

    // Necessary to set up on-line visualization
    m_lattice->copy_data_to_output_buffer();
    m_lattice->post_process();

    // Calculate the number of particles to revert in the context of body force in order to
    // accelerate the flow.
    m_forcing = m_lattice->get_initial_forcing();

    // Setup visualization pipeline
    this->setup_visual();

    // Setup UI
    this->setup_ui();

    m_ui->centralWidget->show();
}

void PipeView::setup_karman()
{
    // TODO
    cout << "Setup Karman vortex street case..." << endl;
}

void PipeView::quit()
{
    qApp->quit();
}

void PipeView::show_about()
{
   system("firefox https://keva92.github.io/lgca/");
}

} // namespace lgca
