/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "writeMeshToFile.h"

#include <vector>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/XDMF/XdmfHdfWriter.h"
#include "MeshLib/Mesh.h"

namespace MeshLib::IO
{
int writeMeshToFile(const MeshLib::Mesh& mesh,
                    std::filesystem::path const& file_path,
                    [[maybe_unused]] std::set<std::string>
                        variable_output_names)
{
    if (file_path.extension().string() == ".msh")
    {
        MeshLib::IO::Legacy::MeshIO meshIO;
        meshIO.setMesh(&mesh);
        BaseLib::IO::writeStringToFile(meshIO.writeToString(), file_path);
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
        std::vector<std::reference_wrapper<const MeshLib::Mesh>> meshes;
        const std::reference_wrapper<const MeshLib::Mesh> mr = mesh;
        meshes.push_back(mr);
        MeshLib::IO::XdmfHdfWriter(std::move(meshes), file_path, 0, 0.0,
                                   variable_output_names, true);
        return 0;
    }
    ERR("writeMeshToFile(): Unknown file extension '{:s}'. Can not write file "
        "'{:s}'.",
        file_path.extension().string(), file_path.string());
    return 0;
}
}  // namespace MeshLib::IO
