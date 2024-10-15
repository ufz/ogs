/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on October 9, 2024, 6:03 PM
 */

#pragma once

namespace MeshLib
{
class Mesh;
}

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;

void checkMPLPhasesForSinglePhaseFlow(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map);
}  // namespace MaterialPropertyLib
