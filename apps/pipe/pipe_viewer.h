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

#ifndef LGCA_PIPE_VIEWER_H_
#define LGCA_PIPE_VIEWER_H_

#include "lgca_common.h"

#include <vtkImageData.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkLookupTable.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#include <QMainWindow>

// Forward Qt class declarations
namespace Ui {
    class PipeView;
}

namespace lgca {

// Forward declarations
template<Model model> class IoVti;
template<Model model> class Lattice;

class PipeView : public QMainWindow
{
    Q_OBJECT

public:

    friend class PipeRunnable;

    explicit PipeView(QWidget *parent = 0);
    ~PipeView();

signals:

public slots:

    void run();
    void stop();

    void rescale();

    void view_cell_density();
    void view_cell_momentum();
    void view_mean_density();
    void view_mean_momentum();

private slots:

    // Collection of private slots to respond to the user activating any of our menu entries
    void setup_single();
    void setup_diffusion();
    void setup_pipe();
    void setup_karman();

    void quit();

    void show_about();

private:

    // Setup visualization pipeline
    void setup_visual();

    // Setup UI
    void setup_ui();

    // Simulation parameters
    static constexpr Model        MODEL         = Model::FHP_III;
    static constexpr unsigned int PP_INTERVAL   = 5;
    static constexpr int          CG_RADIUS     = 10;
           const     string       OUTPUT_DIR    = "./";
           const     string       OUTPUT_FORMAT = "png";

    // Simulation variables
    size_t            m_steps;
    int               m_mnups;
    size_t            m_num_particles;
    std::vector<Real> m_mean_velocity;
    int               m_forcing;
    Real              m_Re = 80.0; // Reynolds number
    Real              m_Ma = 0.3;  // Mach number

    Ui::PipeView*   m_ui;

    Lattice<MODEL>* m_lattice;
    IoVti  <MODEL>* m_vti_io_handler;

    // Data visualization pipeline objects
    vtkImageDataGeometryFilter* m_geom_filter;
    vtkPolyDataMapper*          m_mapper;
    vtkActor*                   m_actor;
    vtkRenderer*                m_ren;
    vtkScalarBarActor*          m_scalar_bar;
    vtkTextProperty*            m_scalar_bar_txt;
    vtkLookupTable*             m_lut;
    vtkRenderWindow*            m_ren_win;
    vtkWindowToImageFilter*     m_png_filter;
    vtkPNGWriter*               m_png_writer;

}; // class PipeView

} // namespace lgca

#endif /* LGCA_PIPE_VIEWER_H_ */
