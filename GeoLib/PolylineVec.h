/**
 * \file PolylineVec.h
 *
 *  Created on 2010-02-09 by Thomas Fischer
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
