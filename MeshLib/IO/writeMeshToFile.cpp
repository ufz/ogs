/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "writeMeshToFile.h"

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/XDMF/Xdmf3Writer.h"
#include "MeshLib/IO/XDMF/transformData.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
namespace IO
{
int writeMeshToFile(const MeshLib::Mesh& mesh,
                    std::filesystem::path const& file_path)
{
    if (file_path.extension().string() == ".msh")
    {
        MeshLib::IO::Legacy::MeshIO meshIO;
        meshIO.setMesh(&mesh);
        meshIO.writeToFile(file_path);
        return 0;
    }
    if (file_path.extension().string() == ".vtu")
    {
        MeshLib::IO::VtuInterface writer(&mesh);
        auto const result = writer.writeToFile(file_path);
        if (!result)
        {
            ERR("writeMeshToFile(): Could not write mesh to '{:s}'.",
                file_path.string());
            return -1;
        }
        return 0;
    }
    if (file_path.extension().string() == ".xdmf")
    {
        if (auto const result = writeXdmf3(mesh, file_path); !result)
        {
            ERR("writeMeshToFile(): Could not write mesh to '{:s}'.",
                file_path.string());
            return -1;
        }
        return 0;
    }
    ERR("writeMeshToFile(): Unknown mesh file format in file {:s}.",
        file_path.string());
    return -1;

}

} // end namespace IO
} // end namespace MeshLib
