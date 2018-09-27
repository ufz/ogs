/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace ProcessLib
{
struct ParameterBase;
class SourceTerm;

std::unique_ptr<SourceTerm> createVolumetricSourceTerm(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& source_term_mesh,
    NumLib::LocalToGlobalIndexMap const& source_term_dof_table,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order, unsigned const shapefunction_order,
    int const variable_id, int const component_id);

}   // namespace ProcessLib
