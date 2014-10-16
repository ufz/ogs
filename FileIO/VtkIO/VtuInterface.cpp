/**
 * \file
 * \author Lars Bilke
 * \date   2014-09-25
 * \brief  Implementation of the BoostVtuInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VtuInterface.h"

#include "logog/include/logog.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>

#include "Mesh.h"
#include "MeshGenerators/VtkMeshConverter.h"

namespace FileIO
{

VtuInterface::VtuInterface() :
	_mesh(nullptr), _use_compressor(false)
{
}

VtuInterface::~VtuInterface()
{}

MeshLib::Mesh* VtuInterface::readVTUFile(const std::string &file_name)
{
	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
		vtkXMLUnstructuredGridReader::New();

	reader->SetFileName(file_name.c_str());
	reader->Update();

	vtkUnstructuredGrid* vtkGrid = reader->GetOutput();

	return MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtkGrid);
}

void VtuInterface::setMesh(const MeshLib::Mesh* mesh)
{
	if (!mesh)
	{
		ERR("VtuInterface::write(): No mesh specified.");
		return;
	}
	this->_mesh = const_cast<MeshLib::Mesh*>(mesh);
};

bool VtuInterface::write()
{
	if(!_mesh)
	{
		ERR("VtuInterface::write(): No mesh specified.");
		return false;
	}

	// TODO write with MappedMesh
	return false;
}

}
