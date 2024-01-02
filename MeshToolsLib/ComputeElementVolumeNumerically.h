/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 26, 2023, 5:52 PM
 */

#pragma once

namespace MeshLib
{
class Element;
}

namespace MeshToolsLib
{

double computeElementVolumeNumerically(MeshLib::Element const& e);
}  // namespace MeshToolsLib
