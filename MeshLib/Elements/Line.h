/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Line class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LINE_H_
#define LINE_H_

#include "TemplateElement.h"
#include "Edge.h"
#include "LineRule2.h"

namespace MeshLib {

typedef TemplateElement<Edge, LineRule2> Line;

}


#endif /* LINE_H_ */
