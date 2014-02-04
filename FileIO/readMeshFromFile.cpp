/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-27
 * \brief  Implementation of readMeshFromFile function.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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

// BaseLib
#include "FileTools.h"
#include "StringTools.h"

// MeshLib
#include "Mesh.h"

// FileIO
#include "Legacy/MeshIO.h"
#include "XmlIO/Boost/BoostVtuInterface.h"
#include "readMeshFromFile.h"

namespace FileIO
{
MeshLib::Mesh* readMeshFromFile(const std::string &file_name)
{
	if (BaseLib::hasFileExtension("msh", file_name))
	{
		FileIO::MeshIO meshIO;
		return meshIO.loadMeshFromFile(file_name);
	}

	if (BaseLib::hasFileExtension("vtu", file_name))
		return FileIO::BoostVtuInterface::readVTUFile(file_name);

	ERR("readMeshFromFile(): Unknown mesh file format in file %s.", file_name.c_str());
	return nullptr;
}
} // end namespace FileIO
