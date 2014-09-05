/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-10-30
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

#include "SystemTools.h"

TEST(BaseLib, EndianLittle) {
    bool isLittle = false;
    int x = 0x00000001;
    if (*(char*)&x)
        isLittle = true;              //am little

	ASSERT_EQ (isLittle, BaseLib::IsLittleEndian());
}

