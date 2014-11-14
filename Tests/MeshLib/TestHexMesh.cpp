/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include "gtest/gtest.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Hex.h"


TEST(MeshLib, HexTests)
{
	MeshLib::Node** nodes = new MeshLib::Node*[20];
	MeshLib::Hex hex(nodes);
	MeshLib::Hex20 hex20(nodes);
}

