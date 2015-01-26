/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Tri class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TRI_H_
#define TRI_H_

#include "TemplateElement.h"
#include "Face.h"
#include "TriRule3.h"

namespace MeshLib {

typedef TemplateElement<Face,TriRule3> Tri;

}

#endif /* TRI_H_ */
