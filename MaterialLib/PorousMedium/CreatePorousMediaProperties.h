/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Parameter.h"
#include "PorousMediaProperties.h"

namespace MeshLib
{
class Mesh;
}

namespace MaterialLib
{
namespace PorousMedium
{
PorousMediaProperties createPorousMediaProperties(
    MeshLib::Mesh& mesh, BaseLib::ConfigTree const& configs,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}
}
