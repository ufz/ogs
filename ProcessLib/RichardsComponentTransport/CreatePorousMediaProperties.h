/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "PorousMediaProperties.h"

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
namespace RichardsComponentTransport
{
PorousMediaProperties createPorousMediaProperties(
    MeshLib::Mesh& mesh, BaseLib::ConfigTree const& porous_medium_configs,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters);
}
}
