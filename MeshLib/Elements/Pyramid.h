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

#pragma once

#include "TemplateElement.h"
#include "PyramidRule5.h"
#include "PyramidRule13.h"

extern template class MeshLib::TemplateElement<MeshLib::PyramidRule13>;
extern template class MeshLib::TemplateElement<MeshLib::PyramidRule5>;

namespace MeshLib {
using Pyramid = TemplateElement<MeshLib::PyramidRule5>;
using Pyramid13 = TemplateElement<MeshLib::PyramidRule13>;
}
