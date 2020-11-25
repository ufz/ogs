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

namespace MeshLib::IO
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
        auto writer =
            std::make_unique<MeshLib::IO::Xdmf3Writer>(MeshLib::IO::Xdmf3Writer(
                file_path, MeshLib::IO::transformGeometry(mesh),
                MeshLib::IO::transformTopology(mesh),
                MeshLib::IO::transformAttributes(mesh), 0));
        // \TODO Errorhandling, Change data model into static and time depended,
        // then is is not neccessary to give time step 0 a special treatment
        // here
        writer->writeStep(0, 0);

        return 0;
    }
    ERR("writeMeshToFile(): Unknown mesh file format in file {:s}.",
        file_path.string());
    return -1;
}

}  // namespace MeshLib::IO
