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

#ifndef LGCA_VIEWER_H_
#define LGCA_VIEWER_H_

#include "lgca_common.h"

#include <QMainWindow>

// Forward Qt class declarations
namespace Ui {
    class LgcaView;
}

namespace lgca {

// Forward declarations
class IoVti;
class Lattice;

class LgcaView : public QMainWindow
{
    Q_OBJECT

public:

    explicit LgcaView(QWidget *parent = 0);
    ~LgcaView();

signals:

public slots:

    void run();

private:

    Ui::LgcaView*   m_ui;
    Lattice*        m_lattice;
    IoVti*    m_vtiIoHandler;
};

} // namespace lgca

#endif /* LGCA_VIEWER_H_ */
