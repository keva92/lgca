/*
 * lgca_io_vti.h
 *
 *  Created on: Oct 19, 2017
 *      Author: Kerstin Vater
 * Description:
 */

#ifndef LGCA_IO_VTI_H_
#define LGCA_IO_VTI_H_

#include "lgca_common.h"

#include <vtkImageData.h>

// Forward declarations
class Lattice;

namespace lgca {

class IoVti
{
public:

    IoVti(Lattice* lattice);
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
