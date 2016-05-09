/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-27
 * \brief  Implementation of readMeshFromFile function.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file readMeshFromFile.cpp
 * @date 2012-09-27
 * @author Karsten Rink
 */

#include "readMeshFromFile.h"

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <boost/algorithm/string/erase.hpp>

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "MeshLib/Mesh.h"

#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

#ifdef USE_PETSC
#include "FileIO/MPI_IO/NodePartitionedMeshReader.h"
#include "MeshLib/NodePartitionedMesh.h"
#endif

namespace MeshLib
{
namespace IO
{
MeshLib::Mesh* readMeshFromFile(const std::string &file_name)
{
#ifdef USE_PETSC
	FileIO::NodePartitionedMeshReader read_pmesh(PETSC_COMM_WORLD);
	const std::string file_name_base = BaseLib::dropFileExtension(file_name);
	return read_pmesh.read(file_name_base);
#else
	if (BaseLib::hasFileExtension("msh", file_name))
	{
		MeshLib::IO::Legacy::MeshIO meshIO;
		return meshIO.loadMeshFromFile(file_name);
	}

	if (BaseLib::hasFileExtension("vtu", file_name))
		return MeshLib::IO::VtuInterface::readVTUFile(file_name);

	ERR("readMeshFromFile(): Unknown mesh file format in file %s.", file_name.c_str());
	return nullptr;
#endif
}
} // end namespace IO
} // end namespace MeshLib
