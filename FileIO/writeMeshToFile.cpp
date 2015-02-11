/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "writeMeshToFile.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "FileTools.h"
#include "StringTools.h"

// MeshLib
#include "Mesh.h"

// FileIO
#include "Legacy/MeshIO.h"
#ifndef OGS_DONT_USE_VTK
#include "FileIO/VtkIO/VtuInterface.h"
#endif

namespace FileIO
{
void writeMeshToFile(const MeshLib::Mesh &mesh, const std::string &file_name)
{
	if (BaseLib::hasFileExtension("msh", file_name))
	{
		Legacy::MeshIO meshIO;
		meshIO.setMesh(&mesh);
		meshIO.writeToFile(file_name);
#ifndef OGS_DONT_USE_VTK
	} else if (BaseLib::hasFileExtension("vtu", file_name)) {
		FileIO::VtuInterface writer(&mesh);
		writer.writeToFile(file_name);
#endif
	} else {
		ERR("writeMeshToFile(): Unknown mesh file format in file %s.", file_name.c_str());
	}
}

} // end namespace FileIO
