/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Tri class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TRI_H_
#define TRI_H_

#include "TemplateElement.h"
#include "TriRule3.h"
#include "TriRule6.h"


extern template class MeshLib::TemplateElement<MeshLib::TriRule3>;
extern template class MeshLib::TemplateElement<MeshLib::TriRule6>;

namespace MeshLib {

typedef TemplateElement<TriRule3> Tri;
typedef TemplateElement<TriRule6> Tri6;

}

#endif /* TRI_H_ */
