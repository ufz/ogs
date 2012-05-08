/**
 * Tet.cpp
 *
 *      Date: 2012/05/03
 *      Author: KR
 */

#include "FemElem.h"

namespace MeshLib {

FemElem::FemElem()
	: _centroid(GEOLIB::Point(0,0,0))
{
}

FemElem::FemElem(const FemElem &elem)
	: _centroid(elem.getCentroid())
{
}

FemElem::~FemElem()
{
}

}

