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
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Hex.h"


TEST(MeshLib, HexTests)
{
	MeshLib::Node** nodes8 = new MeshLib::Node*[8];
	for (int i=0; i<8; i++)
		nodes8[i] = new MeshLib::Node(0, 0, 0);
	MeshLib::Hex hex(nodes8);
	MeshLib::Node** nodes20 = new MeshLib::Node*[20];
	for (int i=0; i<20; i++)
		nodes20[i] = new MeshLib::Node(0, 0, 0);
	MeshLib::Hex20 hex20(nodes20);
}

