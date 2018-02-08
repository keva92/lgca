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

#include "diffusion_viewer.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QStyleFactory>
#include <QVTKOpenGLWidget.h>

// Main function (CPU/host code)
int main(int argc, char **argv) {

    // Needed to ensure appropriate OpenGL context is created for VTK rendering
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());

    // QT Stuff
    QApplication app(argc, argv);
    QApplication::setStyle("fusion");

    QPalette darkPalette;
    darkPalette.setColor(QPalette::Window,          QColor(53,53,53));
    darkPalette.setColor(QPalette::WindowText,      Qt::white);
    darkPalette.setColor(QPalette::Base,            QColor(25,25,25));
    darkPalette.setColor(QPalette::AlternateBase,   QColor(53,53,53));
    darkPalette.setColor(QPalette::ToolTipBase,     Qt::white);
    darkPalette.setColor(QPalette::ToolTipText,     Qt::white);
    darkPalette.setColor(QPalette::Text,            Qt::white);
    darkPalette.setColor(QPalette::Button,          QColor(53,53,53));
    darkPalette.setColor(QPalette::ButtonText,      Qt::white);
    darkPalette.setColor(QPalette::BrightText,      Qt::red);
    darkPalette.setColor(QPalette::Link,            QColor(42, 130, 218));
    darkPalette.setColor(QPalette::Highlight,       QColor(42, 130, 218));
    darkPalette.setColor(QPalette::HighlightedText, Qt::black);

    QApplication::setPalette(darkPalette);
    qApp->setStyleSheet("QToolTip { color: #ffffff; background-color: #2a82da; border: 1px solid white; }");

    lgca::DiffusionView viewer;
    viewer.show();

    return app.exec();
}
