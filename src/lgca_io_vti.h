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
class Lattice;

class IoVti
{
public:

    IoVti(Lattice* lattice, const std::string scalars);
    virtual ~IoVti() { mImageData->Delete(); }

    // Update image data object.
    void update();

    // Write current image data to file.
    void write(const unsigned int step);

          Lattice* lattice()       { assert(mLattice); return mLattice; }
    const Lattice* lattice() const { assert(mLattice); return mLattice; }

          vtkImageData* image()       { assert(mImageData); return mImageData; }
    const vtkImageData* image() const { assert(mImageData); return mImageData; }


private:

    Lattice*        mLattice;
    vtkImageData*   mImageData;

}; // class IoVti

} // namespace lgca

#endif /* LGCA_IO_VTI */
