// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/ConfigTree-fwd.h"
#include "Output.h"

namespace MeshLib
{
class Mesh;
}  // namespace MeshLib

namespace ProcessLib
{
struct SubmeshResiduumOutputConfig
{
    //! Controls the output on the #meshes.
    Output output;

    //! The submeshes for whom additional residuum output should
    //! be produced.
    //!
    //! \attention These meshes must be a non-overlapping cover of the entire
    //! simulation domain (bulk mesh)!
    std::vector<std::reference_wrapper<MeshLib::Mesh>> meshes;
};

SubmeshResiduumOutputConfig createSubmeshResiduumOutputConfig(
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes);
}  // namespace ProcessLib
