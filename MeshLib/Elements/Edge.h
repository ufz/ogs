/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Edge class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EDGE_H_
#define EDGE_H_

#include "TemplateEdge.h"

namespace MeshLib {

typedef TemplateEdge<2,CellType::EDGE2> Edge;

}


#endif /* EDGE_H_ */
