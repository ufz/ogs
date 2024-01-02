/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <string>
#include <vector>

#include "BaseLib/Error.h"

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
class SubmeshAssemblySupport
{
public:
    /// Initializes the assembly on submeshes
    ///
    /// \param meshes the submeshes on whom the assembly shall proceed.
    ///
    /// \attention \c meshes must be a must be a non-overlapping cover of the
    /// entire simulation domain (bulk mesh)!
    ///
    /// \return The names of the residuum vectors that will be assembled.
    virtual std::vector<std::string> initializeAssemblyOnSubmeshes(
        std::vector<std::reference_wrapper<MeshLib::Mesh>> const& meshes)
    {
        DBUG(
            "Default implementation of initializeSubmeshAssembly(). Doing "
            "nothing.");

        if (!meshes.empty())
        {
            OGS_FATAL(
                "Submesh residuum assembly is not implemented for this "
                "process. You can avoid this error message, e.g., by removing "
                "<submesh_residuum_output> from the prj file.");
        }

        return {};
    }

    virtual ~SubmeshAssemblySupport() = default;
};
}  // namespace ProcessLib
