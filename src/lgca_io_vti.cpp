/*
 * lgca_io_vti.cpp
 *
 *  Created on: Oct 19, 2017
 *      Author: Kerstin Vater
 * Description:
 */

#include "lgca_io_vti.h"

#include "lattice.h"

#include "vtkImageImport.h"
#include "vtkImageViewer.h"
#include "vtkRenderer.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"

namespace lgca {

IoVti::IoVti(Lattice* lattice) : mLattice(lattice)
{
//    float cImage[4*4];
//    float value = 0;
//    for(unsigned int row = 0; row < 4; ++row)
//    {
//        for(unsigned int col = 0; col < 4; ++col)
//        {
//            cImage[row * 4 + col] = value;
//            value += 10;
//        }
//    }
//    vtkImageImport* importer = vtkImageImport::New();
//    importer->SetDataSpacing(1, 1, 1);
//    importer->SetDataOrigin(0, 0, 0);
//    importer->SetWholeExtent(0, 4-1, 0, 4-1, 0, 0);
//    importer->SetDataExtentToWholeExtent();
//    importer->SetDataScalarTypeToUnsignedChar();
//    importer->SetNumberOfScalarComponents(1);
//    importer->SetImportVoidPointer(cImage);
//    importer->Update();
//    mImageData = importer->GetOutput();

    mImageData = vtkImageData::New();
    assert(mLattice);
    assert(mImageData);

    mImageData->Initialize();
    mImageData->SetDimensions(mLattice->get_n_x(), mLattice->get_n_y(), 1);
    mImageData->AllocateScalars(/*dataType=*/VTK_FLOAT, /*numComponents=*/1);

//    vtkFloatArray* farray = vtkFloatArray::New();
//    farray->SetNumberOfComponents(1);
//    farray->SetArray((float*)(mLattice->mean_density()), mLattice->get_n_x()*mLattice->get_n_y(), 1);
//    mImageData->GetPointData()->SetScalars(farray);
//    farray->Delete();

    this->update();
}

void
IoVti::update()
{
    assert(mLattice);
    assert(mImageData);

    int* dims = mImageData->GetDimensions();
    assert(dims[0] == mLattice->get_n_x());
    assert(dims[1] == mLattice->get_n_y());
    assert(dims[2] == 1);

    for (int z = 0; z < dims[2]; z++)
        for (int y = 0; y < dims[1]; y++)
            for (int x = 0; x < dims[0]; x++)
            {
                float* pixel = static_cast<float*>(mImageData->GetScalarPointer(x,y,z));
                pixel[0] = mLattice->cell_density(x, y); // TODO Insert reasonable data from mLattice
            }
}

} // namespace lgca
