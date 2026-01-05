// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Surface.h"
#include "TemplateVec.h"

namespace GeoLib
{

/**
 * Class SurfaceVec encapsulate a std::vector of Surfaces
 * and a name.
 * */

using SurfaceVec = TemplateVec<GeoLib::Surface>;

}  // namespace GeoLib
