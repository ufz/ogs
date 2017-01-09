/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TestPermeabilityModel.cpp
 *
 * Created on August 17, 2016, 4:00 PM
 */
#include <gtest/gtest.h>

#include "BaseLib/ConfigTree.h"

#include "Tests/TestTools.h"

#include "MaterialLib/PorousMedium/Permeability/createPermeabilityModel.h"

using namespace MaterialLib;
using namespace MaterialLib::PorousMedium;

using Matrix = Eigen::MatrixXd;

Matrix createTestPermeabilityModel(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("permeability");
    return MaterialLib::PorousMedium::createPermeabilityModel(sub_config);
}

class MaterialPermeabilityTest3D : public ::testing::Test
{
public:
    MaterialPermeabilityTest3D() : _dphi(3, 4)
    {
        // Mock of the gradient of shape functions
        _dphi(0, 0) = 1.0;
        _dphi(0, 1) = 1.0;
        _dphi(0, 2) = 1.0;
        _dphi(0, 3) = 1.0;
        _dphi(1, 0) = 2.0;
        _dphi(1, 1) = 2.0;
        _dphi(1, 2) = 2.0;
        _dphi(1, 3) = 2.0;
        _dphi(2, 0) = 3.0;
        _dphi(2, 1) = 3.0;
        _dphi(2, 2) = 3.0;
        _dphi(2, 3) = 3.0;

        _trans_dphi = _dphi.transpose();
    }

protected:
    Matrix _dphi;
    Matrix _trans_dphi;
};

TEST_F(MaterialPermeabilityTest3D, checkAnisotropicPermeability3D)
{
    const char xml[] =
        "<permeability>"
        "   <values> 2.e-10 0. 0. 0. 3.e-10 0. 0. 0. 4.0e-10 </values> "
        "</permeability>";
    Matrix anisK(createTestPermeabilityModel(xml));

    ASSERT_EQ(3, anisK.rows());

    Matrix laplacian = _trans_dphi * anisK * _dphi;

    ASSERT_NEAR(5.e-9, laplacian(0, 0), 1.e-16);
    ASSERT_NEAR(5.e-9, laplacian(2, 2), 1.e-16);
}

TEST_F(MaterialPermeabilityTest3D, checkIsotropicPermeability3D)
{
    const char xml[] =
        "<permeability>"
        "   <values> 2.e-10</values> "
        "</permeability>";
    Matrix K(createTestPermeabilityModel(xml));

    ASSERT_EQ(1, K.size());

    const double k = K(0, 0);

    Matrix laplacian = _trans_dphi * k * _dphi;

    ASSERT_NEAR(2.8e-9, laplacian(0, 0), 1.e-16);
    ASSERT_NEAR(2.8e-9, laplacian(2, 2), 1.e-16);
}

class MaterialPermeabilityTest2D : public ::testing::Test
{
public:
    MaterialPermeabilityTest2D() : _dphi(2, 3)
    {
        // Mock of the gradient of shape functions
        _dphi(0, 0) = 1.0;
        _dphi(0, 1) = 1.0;
        _dphi(0, 2) = 1.0;
        _dphi(1, 0) = 2.0;
        _dphi(1, 1) = 2.0;
        _dphi(1, 2) = 2.0;

        _trans_dphi = _dphi.transpose();
    }

protected:
    Matrix _dphi;
    Matrix _trans_dphi;
};

TEST_F(MaterialPermeabilityTest2D, checkAnisotropicPermeability2D)
{
    const char xml[] =
        "<permeability>"
        "   <values> 2.e-10 0. 0. 3.e-10</values> "
        "</permeability>";
    Matrix anisK(createTestPermeabilityModel(xml));

    ASSERT_EQ(2, anisK.rows());

    Matrix laplacian = _trans_dphi * anisK * _dphi;

    ASSERT_NEAR(1.4e-9, laplacian(0, 0), 1.e-16);
    ASSERT_NEAR(1.4e-9, laplacian(2, 2), 1.e-16);
}

TEST_F(MaterialPermeabilityTest2D, checkIsotropicPermeability2D)
{
    const char xml[] =
        "<permeability>"
        "   <values> 2.e-10</values> "
        "</permeability>";
    Matrix K(createTestPermeabilityModel(xml));

    ASSERT_EQ(1, K.size());

    const double k = K(0, 0);

    Matrix laplacian = _trans_dphi * k * _dphi;

    ASSERT_NEAR(1.0e-9, laplacian(0, 0), 1.e-16);
    ASSERT_NEAR(1.0e-9, laplacian(2, 2), 1.e-16);
}
