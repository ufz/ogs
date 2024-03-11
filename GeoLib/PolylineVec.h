/**
 * \file
 * \author Thomas Fischer
 * \date   2010-02-09
 * \brief  Definition of the PolylineVec class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
