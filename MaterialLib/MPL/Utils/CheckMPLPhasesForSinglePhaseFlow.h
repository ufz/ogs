// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
