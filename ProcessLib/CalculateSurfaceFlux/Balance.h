/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <logog/include/logog.hpp>

#include <memory>
#include <string>

#include "MeshLib/IO/readMeshFromFile.h"
// TODO used for output, if output classes are ready this has to be changed
#include "MeshLib/IO/writeMeshToFile.h"

namespace ProcessLib
{
struct Balance
{
    Balance(std::string&& balance_mesh_name,
            std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
            std::string&& balance_property_vector_name,
            std::string&& balance_output_mesh_file_name)
        : mesh_name(std::move(balance_mesh_name)),
          property_vector_name(std::move(balance_property_vector_name)),
          output_mesh_file_name(std::move(balance_output_mesh_file_name))
    {
        surface_mesh.reset(MeshLib::IO::readMeshFromFile(mesh_name));

        DBUG(
            "read balance meta data:\n\tbalance mesh:\"%s\"\n\tproperty name: "
            "\"%s\"\n\toutput to: \"%s\"",
            mesh_name.c_str(), property_vector_name.c_str(),
            output_mesh_file_name.c_str());
    }

    void save(double const t) const
    {
        // TODO (TomFischer) output, if output classes are ready this has to be
        // changed
        std::string const fname =
            BaseLib::dropFileExtension(output_mesh_file_name) + "_t_" +
            std::to_string(t) + ".vtu";
        MeshLib::IO::writeMeshToFile(*surface_mesh, fname);
    }

    std::unique_ptr<MeshLib::Mesh> surface_mesh;
    std::string const mesh_name;
    std::string const property_vector_name;
    std::string const output_mesh_file_name;
};
}  // namespace ProcessLib
