/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
