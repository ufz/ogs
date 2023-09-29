/**
 * \file
 * \date   Nov 28, 2017
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <map>
#include <memory>

namespace MeshLib
{
class Mesh;
}

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;

class Medium;

MaterialSpatialDistributionMap createMaterialSpatialDistributionMap(
    std::map<int, std::shared_ptr<Medium>> const& media,
    MeshLib::Mesh const& mesh);
}  // namespace MaterialPropertyLib
