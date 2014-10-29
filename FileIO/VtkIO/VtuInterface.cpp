/**
 * \file
 * \author Lars Bilke
 * \date   2014-09-25
 * \brief  Implementation of the VtuInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VtuInterface.h"

#include "logog/include/logog.hpp"

#include <vtkNew.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "FileTools.h"
#include "InSituLib/VtkMappedMeshSource.h"
#include "Mesh.h"
#include "MeshGenerators/VtkMeshConverter.h"

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
	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
		vtkXMLUnstructuredGridReader::New();

	reader->SetFileName(file_name.c_str());
	reader->Update();

	vtkUnstructuredGrid* vtkGrid = reader->GetOutput();

	std::string const mesh_name (BaseLib::extractBaseNameWithoutExtension(file_name));
	return MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtkGrid, mesh_name);
}

bool VtuInterface::writeToFile(std::string const &file_name)
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
	return (vtuWriter->Write() > 0);
}

}
