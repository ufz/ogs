/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Pyramid class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PYRAMID_H_
#define PYRAMID_H_

#include "TemplateElement.h"
#include "PyramidRule5.h"
#include "PyramidRule13.h"

extern template class MeshLib::TemplateElement<MeshLib::PyramidRule13>;
extern template class MeshLib::TemplateElement<MeshLib::PyramidRule5>;

namespace MeshLib {

typedef TemplateElement<PyramidRule5> Pyramid;
typedef TemplateElement<PyramidRule13> Pyramid13;

}

#endif /* PYRAMID_H_ */
