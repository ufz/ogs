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
	MeshLib::Node** nodes8 = new MeshLib::Node*[8];
	MeshLib::Hex hex(nodes8);
	MeshLib::Node** nodes20 = new MeshLib::Node*[20];
	MeshLib::Hex20 hex20(nodes20);
}

