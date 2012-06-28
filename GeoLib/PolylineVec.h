/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
 * \file PolylineVec.h
 *
 * Created on 2010-02-09 by Thomas Fischer
 */

#ifndef POLYLINEVEC_H_
#define POLYLINEVEC_H_

#include "TemplateVec.h"
#include "Polyline.h"

namespace GeoLib {

/**
 * \ingroup GeoLib
 *
 * \brief class PolylineVec encapsulate a std::vector of Polylines
 * additional one can give the vector of polylines a name
 * */

typedef TemplateVec<Polyline> PolylineVec;

} // end namespace

#endif /* POLYLINEVEC_H_ */
