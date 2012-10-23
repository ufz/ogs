/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file testrunner.cpp
 * Created on 2012-04-29 by Lars Bilke
 *
 */

// ** INCLUDES **
#include "gtest.h"

/// Implementation of the googletest testrunner
int main(int argc, char* argv[])
{
	testing::InitGoogleTest ( &argc, argv );
	return RUN_ALL_TESTS();
}
