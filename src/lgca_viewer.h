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

#ifndef LGCA_VIEW_H_
#define LGCA_VIEW_H_

#include <QMainWindow>

class LgcaView : public QMainWindow
{
    Q_OBJECT

public:

  LgcaView();
  ~LgcaView();

public slots:

  virtual void slotOpenFile();
  virtual void slotExit();

protected:

protected slots:

private:

};

#endif /* LGCA_VIEW_H_ */
