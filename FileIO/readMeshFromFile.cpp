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
	MeshLib::Mesh* mesh (NULL);
	std::string suffix (BaseLib::getSuffixFromPath(file_name));

	if (suffix.compare("msh") == 0 || suffix.compare("MSH") == 0)
	{
		FileIO::MeshIO meshIO;
		mesh = meshIO.loadMeshFromFile(file_name);
	}
	else if (suffix.compare("vtu") == 0 || suffix.compare("VTU") == 0)
		mesh = FileIO::RapidVtuInterface::readVTUFile(file_name);
	else
		std::cout << "Unknown mesh file format" << std::endl;

	return mesh;
}

}
