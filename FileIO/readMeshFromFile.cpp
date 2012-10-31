/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file readMeshFromFile.cpp
 *
 * Created on 2012-09-27 by Karsten Rink
 */

#include "readMeshFromFile.h"
#include "Mesh.h"
#include "RapidXmlIO/RapidVtuInterface.h"
#include "Legacy/MeshIO.h"

// BaseLib
#include "StringTools.h"
#include "FileTools.h"

namespace FileIO {

MeshLib::Mesh* readMeshFromFile(const std::string &file_name)
{
	if (BaseLib::hasFileExtension("msh", file_name))
	{
		FileIO::MeshIO meshIO;
		return meshIO.loadMeshFromFile(file_name);
	}

	if (BaseLib::hasFileExtension("vtu", file_name))
	{
		return FileIO::RapidVtuInterface::readVTUFile(file_name);
	}

	std::cout << "Unknown mesh file format" << std::endl;
	return 0;
}

}
