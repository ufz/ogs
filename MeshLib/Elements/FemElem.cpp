/**
 * \file Tet.cpp
 *
 * Created on 2012-05-03 by Karsten Rink
 */

#include "FemElem.h"

namespace MeshLib {

FemElem::FemElem()
	: _centroid(GeoLib::Point(0,0,0))
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

