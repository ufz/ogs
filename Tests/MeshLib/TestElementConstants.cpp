/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include "gtest/gtest.h"

#include "Elements/Quad.h"
#include "Elements/Tri.h"

using namespace MeshLib;

TEST(MeshLib, ElementConstatns)
{
	ASSERT_EQ(2u, Quad::Dimension);
	ASSERT_EQ(4u, Quad::NAllNodes);
	ASSERT_EQ(4u, Quad::NBaseNodes);
}

