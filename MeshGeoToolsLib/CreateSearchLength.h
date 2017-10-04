/**
 * @brief Functionality to build different search length algorithm objects from
 * given config.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace MeshGeoToolsLib
{
class SearchLength;

std::unique_ptr<MeshGeoToolsLib::SearchLength> createSearchLengthAlgorithm(
    BaseLib::ConfigTree const& external_config, MeshLib::Mesh const& mesh);
}  // end namespace MeshGeoToolsLib
