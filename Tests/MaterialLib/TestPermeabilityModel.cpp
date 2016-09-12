/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TestPermeabilityModel.cpp
 *
 * Created on August 17, 2016, 4:00 PM
 */
#include <gtest/gtest.h>

#include "TestTools.h"

#include "MaterialLib/PorousMedium/Permeability/createPermeabilityModel.h"

namespace
{
using namespace MaterialLib;
using namespace MaterialLib::PorousMedium;

// Mock of the gradient of shape functions
struct MockGradPhi
{
    MockGradPhi()
    {
        dphi.resize(3, 3);
        dphi(0, 0) = 1.0;
        dphi(0, 1) = 1.0;
        dphi(0, 2) = 1.0;
        dphi(1, 0) = 2.0;
        dphi(1, 1) = 2.0;
        dphi(1, 2) = 2.0;
        dphi(2, 0) = 3.0;
        dphi(2, 1) = 3.0;
        dphi(2, 2) = 3.0;

        trans_dphi = dphi.transpose();
    }

    CoefMatrix dphi;
    CoefMatrix trans_dphi;
};

static MockGradPhi mock_grad_phi;

CoefMatrix createTestPermeabilityModel(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("permeability");
    return MaterialLib::PorousMedium::createPermeabilityModel(sub_config);
}

TEST(Material, checkAnisotropicPermeability)
{
    const char xml[] =
        "<permeability>"
        "   <values> 2.e-10 0. 0. 0. 3.e-10 0. 0. 0. 4.0e-10 </values> "
        "</permeability>";
    CoefMatrix anisK(createTestPermeabilityModel(xml));

    ASSERT_EQ(3, anisK.rows());

    CoefMatrix laplacian =
        mock_grad_phi.trans_dphi * anisK * mock_grad_phi.dphi;

    ASSERT_NEAR(5.e-9, laplacian(0, 0), 1.e-16);
    ASSERT_NEAR(5.e-9, laplacian(2, 2), 1.e-16);
}

TEST(Material, checkIsotropicPermeability)
{
    const char xml[] =
        "<permeability>"
        "   <values> 2.e-10</values> "
        "</permeability>";
    CoefMatrix K(createTestPermeabilityModel(xml));

    ASSERT_EQ(1, K.size());

    const double k = K(0, 0);

    CoefMatrix laplacian = mock_grad_phi.trans_dphi * k * mock_grad_phi.dphi;

    ASSERT_NEAR(2.8e-9, laplacian(0, 0), 1.e-16);
    ASSERT_NEAR(2.8e-9, laplacian(2, 2), 1.e-16);
}
}
