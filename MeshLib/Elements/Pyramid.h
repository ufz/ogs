/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Pyramid class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PYRAMID_H_
#define PYRAMID_H_

#include "TemplatePyramid.h"

namespace MeshLib {

typedef TemplatePyramid<5,CellType::PYRAMID5> Pyramid;

}

#endif /* PYRAMID_H_ */
