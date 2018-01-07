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

private:

    // Simulation parameters
    static constexpr Model        MODEL       = Model::FHP;
    static constexpr unsigned int WRITE_STEPS = 20;
    static constexpr int          CG_RADIUS   = 10;             // Coarse graining radius

    // Simulation variables
    int               m_mnups;
    int               m_num_particles;
    std::vector<Real> m_mean_velocity;
    int               m_forcing;
    Real              m_Re = 80.0; // Reynolds number
    Real              m_Ma = 0.2;  // Mach number

    Ui::PipeView*   m_ui;

    Lattice<MODEL>* m_lattice;
    IoVti  <MODEL>* m_vti_io_handler;

    vtkImageDataGeometryFilter* m_geom_filter;
    vtkPolyDataMapper*          m_mapper;
    vtkActor*                   m_actor;
    vtkRenderer*                m_ren;
    vtkScalarBarActor*          m_scalar_bar;
    vtkTextProperty*            m_scalar_bar_txt;
    vtkLookupTable*             m_lut;
    vtkRenderWindow*            m_ren_win;
};

} // namespace lgca

#endif /* LGCA_PIPE_VIEWER_H_ */
