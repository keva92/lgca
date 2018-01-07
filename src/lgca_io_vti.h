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

#ifndef LGCA_IO_VTI_H_
#define LGCA_IO_VTI_H_

#include "lgca_common.h"

#include <vtkImageData.h>

namespace lgca {

// Forward declarations
template<Model model>
class Lattice;

template<Model model_>
class IoVti
{
    using LatticeType = Lattice<model_>;

public:

    IoVti(LatticeType* lattice, const std::string scalars);
    virtual ~IoVti() { m_cell_image_data->Delete(); m_mean_image_data->Delete(); }

    // Set active array for on-line visualization
    void set_scalars(const std::string scalars);

    // Update image data object
    void update();

    // Write current image data to file
    void write(const unsigned int step);

          LatticeType* lattice()       { assert(m_lattice); return m_lattice; }
    const LatticeType* lattice() const { assert(m_lattice); return m_lattice; }

          vtkImageData* cell_image()       { assert(m_cell_image_data); return m_cell_image_data; }
    const vtkImageData* cell_image() const { assert(m_cell_image_data); return m_cell_image_data; }

          vtkImageData* mean_image()       { assert(m_mean_image_data); return m_mean_image_data; }
    const vtkImageData* mean_image() const { assert(m_mean_image_data); return m_mean_image_data; }


private:

    LatticeType*    m_lattice;

    vtkImageData*   m_cell_image_data;
    vtkImageData*   m_mean_image_data;

}; // class IoVti

} // namespace lgca

#endif /* LGCA_IO_VTI */
