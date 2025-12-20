// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ParameterLib
{
struct ParameterBase;
}

namespace ProcessLib
{
class SourceTerm;

std::unique_ptr<SourceTerm> createVolumetricSourceTerm(
    BaseLib::ConfigTree const& config, unsigned const bulk_mesh_dimension,
    MeshLib::Mesh const& source_term_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order, unsigned const shapefunction_order);

}  // namespace ProcessLib
