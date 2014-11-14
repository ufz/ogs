/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-27
 * \brief  Implementation of readMeshFromFile function.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file readMeshFromFile.cpp
 * @date 2012-09-27
 * @author Karsten Rink
 */

// ThirdParty/logog
#include "logog/include/logog.hpp"
#ifdef USE_PETSC
#include <petscksp.h>
#endif

// BaseLib
#include "FileTools.h"
#include "StringTools.h"

// MeshLib
#include "Mesh.h"

// FileIO
#include "Legacy/MeshIO.h"
#include "FileIO/VtkIO/VtuInterface.h"
#include "readMeshFromFile.h"
//
#ifdef USE_PETSC
#include "MPI_IO/NodePartitionedMeshReader.h"
#include "MeshLib/NodePartitionedMesh.h"
#endif


namespace FileIO
{
MeshLib::Mesh* readMeshFromFile(const std::string &file_name)
{
#ifdef USE_PETSC
	NodePartitionedMeshReader read_pmesh;
	MeshLib::NodePartitionedMesh *mesh 
		= read_pmesh.read(PETSC_COMM_WORLD, BaseLib::extractBaseName(file_name));
	if(mesh)
		return mesh;
#else
	if (BaseLib::hasFileExtension("msh", file_name))
	{
		Legacy::MeshIO meshIO;
		return meshIO.loadMeshFromFile(file_name);
	}

	if (BaseLib::hasFileExtension("vtu", file_name))
		return VtuInterface::readVTUFile(file_name);
#endif

	ERR("readMeshFromFile(): Unknown mesh file format in file %s.", file_name.c_str());
	return nullptr;
}
} // end namespace FileIO
