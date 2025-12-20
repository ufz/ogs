// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
