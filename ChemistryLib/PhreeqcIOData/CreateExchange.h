/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <optional>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct ExchangeSite;

std::vector<ExchangeSite> createExchange(
    std::optional<BaseLib::ConfigTree> const& config, MeshLib::Mesh& mesh);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
