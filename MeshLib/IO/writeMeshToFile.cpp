// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
                        output_variable_names,
                    bool const use_compression,
                    int const data_mode)
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
        if (output_variable_names.empty())
        {
            for (auto const& name :
                 mesh.getProperties().getPropertyVectorNames())
            {
                output_variable_names.insert(name);
            }
        }
        MeshLib::IO::VtuInterface writer(&mesh, output_variable_names,
                                         data_mode, use_compression);
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
                                   output_variable_names, use_compression, 1,
                                   1048576);
        return 0;
    }
    ERR("writeMeshToFile(): Unknown file extension '{:s}'. Can not write file "
        "'{:s}'.",
        file_path.extension().string(), file_path.string());
    return -1;
}
}  // namespace MeshLib::IO
