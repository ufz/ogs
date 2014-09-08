/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Tet class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TET_H_
#define TET_H_

#include "TemplateTet.h"

namespace MeshLib {

typedef TemplateTet<4,CellType::TET4> Tet;

}

#endif /* TET_H_ */
