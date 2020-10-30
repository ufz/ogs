/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "writeMeshToFile.h"

#include "BaseLib/Logging.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/XDMF/writeXdmf.h"


namespace MeshLib
{
namespace IO
{
int writeMeshToFile(const MeshLib::Mesh &mesh, const std::string &file_name)
{
    if (BaseLib::hasFileExtension(".msh", file_name))
    {
        MeshLib::IO::Legacy::MeshIO meshIO;
        meshIO.setMesh(&mesh);
        meshIO.writeToFile(file_name);
        return 0;
    }
    if (BaseLib::hasFileExtension(".vtu", file_name))
    {
        MeshLib::IO::VtuInterface writer(&mesh);
        auto const result = writer.writeToFile(file_name);
        if (!result)
        {
            ERR("writeMeshToFile(): Could not write mesh to '{:s}'.",
                file_name);
            return -1;
        }
        return 0;
    }
    if (BaseLib::hasFileExtension(".xdmf", file_name))
    {
        if (auto const result = writeXdmf3(mesh, file_name); !result)
        {
            ERR("writeMeshToFile(): Could not write mesh to '{:s}'.",
                file_name);
            return -1;
        }
        return 0;
    }
    ERR("writeMeshToFile(): Unknown mesh file format in file {:s}.", file_name);
    return -1;

}

} // end namespace IO
} // end namespace MeshLib
