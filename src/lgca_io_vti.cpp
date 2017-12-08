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
#include "vtkSOADataArrayTemplate.h"
#include "vtkCellData.h"
#include "vtkXMLImageDataWriter.h"

#include <sstream>

namespace lgca {

IoVti::IoVti(Lattice* lattice) : mLattice(lattice)
{
    mImageData = vtkImageData::New();
    assert(mLattice);
    assert(mImageData);

    mImageData->Initialize();
    mImageData->SetDimensions(mLattice->get_n_x() + 1, mLattice->get_n_y() + 1, 1); // Number of points in each direction

    // Pass pointer to cell density array of the lattice to the image data object
    vtkFloatArray* cell_density = vtkFloatArray::New();
    cell_density->SetName("Cell density");
    cell_density->SetNumberOfComponents(1);
    cell_density->SetArray((float*)(mLattice->cell_density()), mLattice->get_n_x()*mLattice->get_n_y(), /*save=*/1);
    mImageData->GetCellData()->AddArray(cell_density);
    cell_density->Delete();

    // Pass pointer to mean density array of the lattice to the image data object
    vtkFloatArray* mean_density = vtkFloatArray::New();
    mean_density->SetName("Mean density");
    mean_density->SetNumberOfComponents(1);
    mean_density->SetArray((float*)(mLattice->mean_density()), mLattice->get_n_x()*mLattice->get_n_y(), /*save=*/1);
    mImageData->GetCellData()->AddArray(mean_density);
    mean_density->Delete();

    // Pass pointer to cell momentum array of the lattice to the image data object
    vtkSOADataArrayTemplate<float>* cell_momentum = vtkSOADataArrayTemplate<float>::New();
    cell_momentum->SetName("Cell momentum");
    cell_momentum->SetNumberOfComponents(2);
    cell_momentum->SetArray(/*comp=*/0, (float*)(mLattice->cell_momentum()), mLattice->get_n_x()*mLattice->get_n_y(), /*updateMaxId=*/0, /*save=*/1);
    cell_momentum->SetArray(/*comp=*/1, (float*)(mLattice->cell_momentum()) + mLattice->get_n_x()*mLattice->get_n_y(), mLattice->get_n_x()*mLattice->get_n_y(), /*updateMaxId=*/0, /*save=*/1);
    mImageData->GetCellData()->AddArray(cell_momentum);
    cell_momentum->Delete();

    // Set active array for on-line visualization
    mImageData->GetCellData()->SetActiveScalars("Mean density");

    // Mark image data object as modified
    this->update();
}

void
IoVti::update()
{
    assert(mImageData);
    mImageData->Modified();
}

void
IoVti::write(const unsigned int step)
{
    // TODO Set flags which data should be written to file.
    bool write_cell_density  = true;
    bool write_mean_density  = true;
    bool write_cell_momentum = false;
    bool write_mean_momentum = false;
    bool write_cell_velocity = false;
    bool write_mean_velocity = false;

    // Compose file name.
    std::ostringstream filename;
    filename << "../res/res_" << step << ".vti";

    printf("Writing results to file %s...\n", filename.str().c_str());

    vtkXMLImageDataWriter* writer = vtkXMLImageDataWriter::New();
    writer->SetFileName(filename.str().c_str());
    writer->SetInputData(mImageData);
    writer->SetCompressorTypeToLZ4();
    writer->SetDataModeToBinary();
    writer->Write();
    writer->Delete();

    printf("...done.\n");
}

} // namespace lgca
