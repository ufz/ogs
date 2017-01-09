/**
 * \file
 * \author Thomas Fischer
 * \date   2010-02-09
 * \brief  Definition of the SurfaceVec class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TemplateVec.h"
#include "Surface.h"

namespace GeoLib {

/**
 * Class SurfaceVec encapsulate a std::vector of Surfaces
 * and a name.
 * */

typedef TemplateVec<Surface> SurfaceVec;

} // end namespace
