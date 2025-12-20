// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
