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
#include "vtkPointData.h"
#include "vtkXMLImageDataWriter.h"

#include <sstream>

namespace lgca {

template<int num_dir_>
IoVti<num_dir_>::IoVti(LatticeType* lattice, const std::string scalars) : m_lattice(lattice)
{
    assert(m_lattice);

    m_cell_image_data = vtkImageData::New();
    m_mean_image_data = vtkImageData::New();
    assert(m_cell_image_data);
    assert(m_mean_image_data);

    m_cell_image_data->Initialize();
    m_mean_image_data->Initialize();

    m_mean_image_data->SetDimensions(m_lattice->coarse_dim_x(), m_lattice->coarse_dim_y(), 1);         // TODO Use vtkCellDataToPointData?
    m_cell_image_data->SetDimensions(m_lattice->       dim_x() + 1, m_lattice->       dim_y() + 1, 1); // Number of points in each direction

    // Pass pointer to cell density array of the lattice to the image data object
    vtkFloatArray* cell_density = vtkFloatArray::New();
    cell_density->SetName("Cell density");
    cell_density->SetNumberOfComponents(1);
    cell_density->SetArray((float*)(m_lattice->cell_density()), m_lattice->num_cells(), /*save=*/1);
    m_cell_image_data->GetCellData()->AddArray(cell_density);
    cell_density->Delete();

    // Pass pointer to mean density array of the lattice to the image data object
    vtkFloatArray* mean_density = vtkFloatArray::New();
    mean_density->SetName("Mean density");
    mean_density->SetNumberOfComponents(1);
    mean_density->SetArray((float*)(m_lattice->mean_density()), m_lattice->num_coarse_cells(), /*save=*/1);
    m_mean_image_data->GetPointData()->AddArray(mean_density);
    mean_density->Delete();

    // Pass pointer to cell momentum array of the lattice to the image data object
    vtkAOSDataArrayTemplate<float>* cell_momentum = vtkAOSDataArrayTemplate<float>::New();
    cell_momentum->SetName("Cell momentum");
    cell_momentum->SetNumberOfComponents(2);
    cell_momentum->SetArray((float*)(m_lattice->cell_momentum()), /*size=*/2 * m_lattice->num_cells(), /*save=*/1);
    m_cell_image_data->GetCellData()->AddArray(cell_momentum);
    cell_momentum->Delete();

    // Pass pointer to mean momentum array of the lattice to the image data object
    vtkAOSDataArrayTemplate<float>* mean_momentum = vtkAOSDataArrayTemplate<float>::New();
    mean_momentum->SetName("Mean momentum");
    mean_momentum->SetNumberOfComponents(2);
    mean_momentum->SetArray((float*)(m_lattice->mean_momentum()), /*size=*/2 * m_lattice->num_coarse_cells(), /*save=*/1);
    m_mean_image_data->GetPointData()->AddArray(mean_momentum);
    mean_momentum->Delete();

    // Set active array for on-line visualization
    m_cell_image_data->GetCellData()->SetActiveScalars(scalars.c_str());
    m_mean_image_data->GetPointData()->SetActiveScalars(scalars.c_str());

    // Mark image data object as modified
    this->update();
}

template<int num_dir_>
void IoVti<num_dir_>::update()
{
    assert(m_cell_image_data);
    assert(m_mean_image_data);

    m_cell_image_data->Modified();
    m_mean_image_data->Modified();
}

template<int num_dir_>
void IoVti<num_dir_>::write(const unsigned int step)
{
    // TODO Set flags which data should be written to file.
    bool write_cell_density  = true;
    bool write_mean_density  = true;
    bool write_cell_momentum = false;
    bool write_mean_momentum = false;
    bool write_cell_velocity = false;
    bool write_mean_velocity = false;

    // Compose file name
    std::ostringstream cell_filename, mean_filename;
    cell_filename << "../res/cell_res_" << step << ".vti";
    mean_filename << "../res/mean_res_" << step << ".vti";

    printf("Writing results to files %s and %s...\n",
           cell_filename.str().c_str(), mean_filename.str().c_str());

    vtkXMLImageDataWriter* cell_data_writer = vtkXMLImageDataWriter::New();
    cell_data_writer->SetFileName(cell_filename.str().c_str());
    cell_data_writer->SetInputData(m_cell_image_data);
    cell_data_writer->SetCompressorTypeToLZ4();
    cell_data_writer->SetDataModeToBinary();
    cell_data_writer->Write();
    cell_data_writer->Delete();

    vtkXMLImageDataWriter* mean_data_writer = vtkXMLImageDataWriter::New();
    mean_data_writer->SetFileName(mean_filename.str().c_str());
    mean_data_writer->SetInputData(m_mean_image_data);
    mean_data_writer->SetCompressorTypeToLZ4();
    mean_data_writer->SetDataModeToBinary();
    mean_data_writer->Write();
    mean_data_writer->Delete();

    printf("...done.\n");
}

// Explicit instantiations
template class IoVti<4>;
template class IoVti<6>;

} // namespace lgca
