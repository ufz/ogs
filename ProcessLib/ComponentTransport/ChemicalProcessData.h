/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ProcessLib
{
namespace ComponentTransport
{
struct ChemicalProcessData
{
    ChemicalProcessData(
        std::vector<std::vector<GlobalIndexType>>& chemical_system_index_map_)
        : chemical_system_index_map(chemical_system_index_map_)
    {
    }

    std::vector<std::vector<GlobalIndexType>>& chemical_system_index_map;
};
}  // namespace ComponentTransport
}  // namespace ProcessLib
