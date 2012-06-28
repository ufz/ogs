/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SurfaceVec.h
 *
 * Created on 2010-02-09 by Thomas Fischer
 */


#ifndef SURFACEVEC_H_
#define SURFACEVEC_H_

#include "TemplateVec.h"
#include "Surface.h"

namespace GeoLib {

/**
 * Class SurfaceVec encapsulate a std::vector of Surfaces
 * and a name.
 * */

typedef TemplateVec<Surface> SurfaceVec;

} // end namespace

#endif /* SURFACEVEC_H_ */
