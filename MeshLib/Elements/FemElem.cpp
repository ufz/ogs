/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-03
 * \brief  Implementation of the FemElem class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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

