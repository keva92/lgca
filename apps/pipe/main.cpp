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

#include "pipe_viewer.h"

#include "lgca_bitset.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLWidget.h>

// Main function (CPU/host code)
int main(int argc, char **argv) {

    // Needed to ensure appropriate OpenGL context is created for VTK rendering
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());

    // QT Stuff
    QApplication app(argc, argv);
    QApplication::setStyle("fusion");

    lgca::PipeView viewer;
    viewer.show();

    return app.exec();
}
