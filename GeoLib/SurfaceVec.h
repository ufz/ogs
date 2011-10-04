/*
 * \file SurfaceVec.h
 *
 *  Created on: Feb 9, 2010
 *      Author: fischeth
 */


#ifndef SURFACEVEC_H_
#define SURFACEVEC_H_

#include "TemplateVec.h"
#include "Surface.h"

namespace GEOLIB {

/**
 * Class SurfaceVec encapsulate a std::vector of Surfaces
 * and a name.
 * */

typedef TemplateVec<Surface> SurfaceVec;

} // end namespace

#endif /* SURFACEVEC_H_ */
