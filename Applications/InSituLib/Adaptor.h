/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"

namespace MeshLib {
    class Mesh;
}

namespace InSituLib {
void Initialize(BaseLib::ConfigTree const& scripts_config,
                std::string const& path);
void Finalize();
void CoProcess(MeshLib::Mesh const& mesh, double const time,
               unsigned int const timeStep, bool const lastTimeStep);
}
