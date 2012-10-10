/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file Edge.h
 *
 *  Created on  Sep 27, 2012 by Thomas Fischer
 */

#ifndef EDGE_H_
#define EDGE_H_

#include "TemplateEdge.h"

namespace MeshLib {

typedef TemplateEdge<2,FEMElemType::EDGE2> Edge;

}


#endif /* EDGE_H_ */
