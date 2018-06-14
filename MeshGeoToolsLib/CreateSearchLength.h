/**
 * @brief Functionality to build different search length algorithm objects from
 * given config.
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

/// Creates a search length algorithm from the given config.
///
/// In case, that there is no tag for the search length algorithm found in the
/// config, the default SearchLength algorithm is returned.
std::unique_ptr<MeshGeoToolsLib::SearchLength> createSearchLengthAlgorithm(
    BaseLib::ConfigTree const& external_config, MeshLib::Mesh const& mesh);
}  // end namespace MeshGeoToolsLib
