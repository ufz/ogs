/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 22, 2021, 11:52 AM
 */

#pragma once

#include <vector>

namespace MeshLib
{
class Mesh;

bool is2DMeshOnRotatedVerticalPlane(Mesh const& mesh);

};  // namespace MeshLib
