/**
 * \file
 * \author Lars Bilke
 * \date   2014-09-25
 * \brief  Implementation of the VtuInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VtuInterface.h"

#include <vtkNew.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <logog/include/logog.hpp>

#include <boost/algorithm/string/erase.hpp>

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include "BaseLib/FileTools.h"
#include "InSituLib/VtkMappedMeshSource.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/VtkMeshConverter.h"

namespace FileIO
{

VtuInterface::VtuInterface(const MeshLib::Mesh* mesh, int dataMode, bool compress) :
	_mesh(mesh), _data_mode(dataMode), _use_compressor(compress)
{
	if(_data_mode == vtkXMLWriter::Appended)
		ERR("Appended data mode is currently not supported!");
	if(_data_mode == vtkXMLWriter::Ascii && compress)
		WARN("Ascii data cannot be compressed, ignoring compression flag.")
}

VtuInterface::~VtuInterface()
{}

MeshLib::Mesh* VtuInterface::readVTUFile(std::string const &file_name)
{
	if (!BaseLib::IsFileExisting(file_name))
		return nullptr;

	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
		vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName(file_name.c_str());
	reader->Update();

	vtkUnstructuredGrid* vtkGrid = reader->GetOutput();
	if (vtkGrid->GetNumberOfPoints() == 0)
		return nullptr;

	std::string const mesh_name (BaseLib::extractBaseNameWithoutExtension(file_name));
	return MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtkGrid, mesh_name);
}

bool VtuInterface::writeToFile(std::string const &file_name)
{
#ifdef USE_PETSC
	// Also for other approach with DDC.
	// In such case, a MPI_Comm argument is need to this member,
	// and PETSC_COMM_WORLD shoud be replaced with the argument.  
	int mpi_rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
	const std::string file_name_base = boost::erase_last_copy(file_name, ".vtu");

	const std::string file_name_rank = file_name_base + "_"
	                                   + std::to_string(mpi_rank) + ".vtu";
	const bool vtu_status_i = writeVTU(file_name_rank);
	bool vtu_status = false;
    MPI_Allreduce(&vtu_status_i, &vtu_status, 1, MPI_C_BOOL, MPI_LAND, PETSC_COMM_WORLD);

	int mpi_size;
	MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
	bool pvtu_status = false;
	if (mpi_rank == 0)
	{
		pvtu_status = writeVTU(file_name_base + ".pvtu", mpi_size);
	}
	MPI_Bcast(&pvtu_status, 1, MPI_C_BOOL, 0, PETSC_COMM_WORLD);

	return vtu_status && pvtu_status;

#else
	return writeVTU(file_name);
#endif
}

bool VtuInterface::writeVTU(std::string const &file_name, const int num_partitions)
{
	if(!_mesh)
	{
		ERR("VtuInterface::write(): No mesh specified.");
		return false;
	}

	// See http://www.paraview.org/Bug/view.php?id=13382
	if(_data_mode == vtkXMLWriter::Appended)
		WARN("Appended data mode is currently not supported, written file is not valid!");

	vtkNew<InSituLib::VtkMappedMeshSource> vtkSource;
	vtkSource->SetMesh(_mesh);

	vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriter =
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

	vtuWriter->SetInputConnection(vtkSource->GetOutputPort());
	if(_use_compressor)
		vtuWriter->SetCompressorTypeToZLib();

	vtuWriter->SetDataMode(_data_mode);
	if (_data_mode == vtkXMLWriter::Appended)
		vtuWriter->SetEncodeAppendedData(1);
	if (_data_mode == vtkXMLWriter::Ascii)
	{
		// Mapped data structures for OGS to VTK mesh conversion are not fully
		// implemented and doing so is not trivial. Therefore for ascii output
		// the mapped unstructured grid is copied to a regular VTK grid.
		// See http://www.vtk.org/pipermail/vtkusers/2014-October/089400.html
		vtkSource->Update();
		vtkSmartPointer<vtkUnstructuredGrid> tempGrid =
			vtkSmartPointer<vtkUnstructuredGrid>::New();
		tempGrid->DeepCopy(vtkSource->GetOutput());
		vtuWriter->SetInputDataObject(tempGrid);
	}

	vtuWriter->SetFileName(file_name.c_str());
	if (num_partitions > 0)
		vtuWriter->SetNumberOfPieces(num_partitions);

	return (vtuWriter->Write() > 0);
}

}
