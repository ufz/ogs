// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/ConfigTree-fwd.h"

namespace MeshLib {
    class Mesh;
}

namespace InSituLib {
void Initialize(BaseLib::ConfigTree const& scripts_config,
                std::string const& path);
void Finalize();
void CoProcess(MeshLib::Mesh const& mesh, double const time,
               unsigned int const timeStep, bool const lastTimeStep,
               std::string output_directory);
}
