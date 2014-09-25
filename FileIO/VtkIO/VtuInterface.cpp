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

#include "Mesh.h"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>

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

	return nullptr;
}

}
