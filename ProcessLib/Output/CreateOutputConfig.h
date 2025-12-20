// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/ConfigTree-fwd.h"
#include "OutputConfig.h"

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
OutputConfig createOutputConfig(
    const BaseLib::ConfigTree& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes);
}
