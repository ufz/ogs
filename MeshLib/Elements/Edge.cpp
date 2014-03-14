/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Edge class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Edge.h"
#include "Node.h"

#include "MathTools.h"


namespace MeshLib
{
Edge::Edge(unsigned value, unsigned id)
    : Element(value, id), _length(-1.0) // init with invalid value to detect errors
{
}


}

