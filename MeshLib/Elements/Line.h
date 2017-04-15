/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Definition of the Line class.
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
#include "LineRule2.h"
#include "LineRule3.h"

extern template class MeshLib::TemplateElement<MeshLib::LineRule2>;
extern template class MeshLib::TemplateElement<MeshLib::LineRule3>;

namespace MeshLib {
using Line = TemplateElement<MeshLib::LineRule2>;
using Line3 = TemplateElement<MeshLib::LineRule3>;
}
