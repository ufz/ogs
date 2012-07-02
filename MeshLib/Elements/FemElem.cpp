/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElem.cpp
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

