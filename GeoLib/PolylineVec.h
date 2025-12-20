// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Polyline.h"
#include "TemplateVec.h"

namespace GeoLib
{

/**
 * \brief class PolylineVec encapsulate a std::vector of Polylines
 * additional one can give the vector of polylines a name
 * */
using PolylineVec = TemplateVec<GeoLib::Polyline>;

}  // namespace GeoLib
