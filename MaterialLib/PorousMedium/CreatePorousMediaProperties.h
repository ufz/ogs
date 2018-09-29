/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "PorousMediaProperties.h"
#include "BaseLib/ConfigTree.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace MeshLib
{
class Mesh;
}

namespace MaterialLib
{
namespace PorousMedium
{
PorousMediaProperties createPorousMediaProperties(
    MeshLib::Mesh& mesh, BaseLib::ConfigTree const& porous_media_config,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters);
}
}
