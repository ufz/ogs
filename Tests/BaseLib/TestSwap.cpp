/**
 * \file
 * \author Lars Bilke
 * \date   2010-04-29
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "gtest/gtest.h"

TEST(BaseLib, SwapInt) {
	int arg0 = 5;
	int arg1 = 10;
	std::swap(arg0, arg1);
	ASSERT_EQ ( arg0, 10 );
	ASSERT_EQ ( arg1, 5 );
}

TEST(BaseLib, SwapDouble) {
	double arg0 = 5.0;
	double arg1 = 10.0;
	std::swap(arg0, arg1);
	ASSERT_EQ ( arg0, 10.0 );
	ASSERT_EQ ( arg1, 5.0 );
}

TEST(BaseLib, SwapString) {
	std::string arg0 = "5";
	std::string arg1 = "10";
	std::swap(arg0, arg1);
	ASSERT_EQ ( arg0, std::string("10") );
	ASSERT_EQ ( arg1, std::string("5") );
}
