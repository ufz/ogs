/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-27
 * \brief  Implementation of readMeshFromFile function.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file readMeshFromFile.cpp
 * @date 2012-09-27
 * @author Karsten Rink
 */

#ifdef USE_PETSC
#include <petsc.h>
#endif

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "FileTools.h"
#include "StringTools.h"

// MeshLib
#include "Mesh.h"

// FileIO
#include "Legacy/MeshIO.h"
#include "FileIO/VtkIO/VtuInterface.h"
#include "readMeshFromFile.h"
// FileIO : for reading partitioned mesh.
#ifdef USE_PETSC
#include "FileIO/MPI_IO/NodePartitionedMeshReader.h"
#include "MeshLib/NodePartitionedMesh.h"
#endif

namespace FileIO
{
MeshLib::Mesh* readMeshFromFile(const std::string &file_name)
{
#ifdef USE_PETSC
	NodePartitionedMeshReader read_pmesh(PETSC_COMM_WORLD);
	return read_pmesh.read(file_name);
#else
	if (BaseLib::hasFileExtension("msh", file_name))
	{
		Legacy::MeshIO meshIO;
		return meshIO.loadMeshFromFile(file_name);
	}

	if (BaseLib::hasFileExtension("vtu", file_name))
		return VtuInterface::readVTUFile(file_name);

	ERR("readMeshFromFile(): Unknown mesh file format in file %s.", file_name.c_str());
	return nullptr;
#endif
}
} // end namespace FileIO
