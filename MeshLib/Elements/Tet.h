/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Tet class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TET_H_
#define TET_H_

#include "TemplateElement.h"
#include "TetRule4.h"
#include "TetRule10.h"

namespace MeshLib {

typedef TemplateElement<TetRule4> Tet;
typedef TemplateElement<TetRule10> Tet10;

}

#endif /* TET_H_ */
