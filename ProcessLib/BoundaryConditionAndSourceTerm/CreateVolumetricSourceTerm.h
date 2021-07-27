/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
