/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Quad class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QUAD_H_
#define QUAD_H_

#include "TemplateElement.h"
#include "Face.h"
#include "QuadRule4.h"
#include "QuadRule8.h"
#include "QuadRule9.h"

namespace MeshLib
{

typedef TemplateElement<Face, QuadRule4> Quad;
typedef TemplateElement<Face, QuadRule8> Quad8;
typedef TemplateElement<Face, QuadRule9> Quad9;

}

#endif /* QUAD_H_ */
