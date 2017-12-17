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

IoVti::IoVti(Lattice* lattice, const std::string scalars) : mLattice(lattice)
{
    mImageData = vtkImageData::New();
    assert(mLattice);
    assert(mImageData);

    mImageData->Initialize();
    mImageData->SetDimensions(mLattice->dim_x() + 1, mLattice->dim_y() + 1, 1); // Number of points in each direction

    // Pass pointer to cell density array of the lattice to the image data object
    vtkFloatArray* cell_density = vtkFloatArray::New();
    cell_density->SetName("Cell density");
    cell_density->SetNumberOfComponents(1);
    cell_density->SetArray((float*)(mLattice->cell_density()), mLattice->num_cells(), /*save=*/1);
    mImageData->GetCellData()->AddArray(cell_density);
    cell_density->Delete();

    // Pass pointer to mean density array of the lattice to the image data object
    vtkFloatArray* mean_density = vtkFloatArray::New();
    mean_density->SetName("Mean density");
    mean_density->SetNumberOfComponents(1);
    mean_density->SetArray((float*)(mLattice->mean_density()), mLattice->num_cells(), /*save=*/1);
    mImageData->GetCellData()->AddArray(mean_density);
    mean_density->Delete();

    // Pass pointer to cell momentum array of the lattice to the image data object
    vtkSOADataArrayTemplate<float>* cell_momentum = vtkSOADataArrayTemplate<float>::New();
    cell_momentum->SetName("Cell momentum");
    cell_momentum->SetNumberOfComponents(2);
    cell_momentum->SetArray(/*comp=*/0, (float*)(mLattice->cell_momentum()), mLattice->num_cells(), /*updateMaxId=*/0, /*save=*/1);
    cell_momentum->SetArray(/*comp=*/1, (float*)(mLattice->cell_momentum()) + mLattice->num_cells(), mLattice->num_cells(), /*updateMaxId=*/0, /*save=*/1);
    mImageData->GetCellData()->AddArray(cell_momentum);
    cell_momentum->Delete();

    // Set active array for on-line visualization
    mImageData->GetCellData()->SetActiveScalars(scalars.c_str());

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
