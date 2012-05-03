/**
 * Tet.cpp
 *
 *      Date: 2012/05/03
 *      Author: KR
 */

#include "FemElem.h"

namespace MeshLib {

FemElem::FemElem()
	: _centre_of_gravity(GEOLIB::Point(0,0,0))
{
}

FemElem::FemElem(const FemElem &elem)
	: _centre_of_gravity(elem.getCentreOfGravity())
{
}

FemElem::~FemElem()
{
}

}

