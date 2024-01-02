/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifdef OGS_USE_MFRONT

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include "MaterialLib/SolidModels/MFront/MFrontGeneric.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;
namespace MSM = MaterialLib::Solids::MFront;

TEST(MaterialLib_OgsToMFrontConversion, Tensor3D)
{
    MPL::Tensor<3> ogs_tensor;
    ogs_tensor << 11, 12, 13, 21, 22, 23, 31, 32, 33;

    EXPECT_THAT(MSM::ogsTensorToMFrontTensor<3>(ogs_tensor).eval(),
                testing::Pointwise(testing::DoubleEq(),
                                   {11, 22, 33, 12, 21, 13, 31, 23, 32}));
}

TEST(MaterialLib_OgsToMFrontConversion, Tensor2D)
{
    MPL::Tensor<2> ogs_tensor;
    ogs_tensor << 11, 12, 21, 22, 33;

    EXPECT_THAT(MSM::ogsTensorToMFrontTensor<2>(ogs_tensor).eval(),
                testing::Pointwise(testing::DoubleEq(), {11, 22, 33, 12, 21}));
}

#endif
