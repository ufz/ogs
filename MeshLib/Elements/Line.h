/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Line class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LINE_H_
#define LINE_H_

#include "TemplateLine.h"

namespace MeshLib {

typedef TemplateLine<2,CellType::LINE2> Line;
typedef TemplateLine<2,CellType::LINE3> Line3;

}


#endif /* LINE_H_ */
