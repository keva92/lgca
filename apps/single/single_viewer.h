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

#ifndef LGCA_SINGLE_VIEWER_H_
#define LGCA_SINGLE_VIEWER_H_

#include "lgca_common.h"

#include <QMainWindow>

// Forward Qt class declarations
namespace Ui {
    class SingleView;
}

namespace lgca {

// Forward declarations
class IoVti;
class Lattice;

class SingleView : public QMainWindow
{
    Q_OBJECT

public:

    explicit SingleView(QWidget *parent = 0);
    ~SingleView();

signals:



public slots:

    void run();
    void stop();

private:

    Ui::SingleView* m_ui;
    Lattice*        m_lattice;
    IoVti*          m_vti_io_handler;

    int             m_mnups;
    int             m_num_particles;

    const int       m_write_steps = 1; // Number of steps after which results are post-processed and visualized or written to file
};

} // namespace lgca

#endif /* LGCA_SINGLE_VIEWER_H_ */
