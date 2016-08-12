/*!
   \file  TestTensorMaterialParameter.cpp
   \brief Test classes for tensor material parameters.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>

#include <gtest/gtest.h>

#include "MaterialLib/TensorParameter.h"
#include "MaterialLib/Fluid/Permeability/PermeabilityType.h"
#include "MaterialLib/Fluid/Permeability/IntrinsicPermeability.h"

namespace
{
using namespace MaterialLib;
using namespace MaterialLib::Fluid;

using Matrix = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;

// Mock of the gradient of shape functions
struct MockGrabPhi
{
    MockGrabPhi()
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

    Matrix dphi;
    Matrix trans_dphi;
};

static MockGrabPhi mock_grad_phi;

TEST(Material, checkConstantTensor)
{
    Matrix anis_k(3, 3);
    anis_k.setZero();
    anis_k(0, 0) = 2.e-10;
    anis_k(1, 1) = 3.e-10;
    anis_k(2, 2) = 4.e-10;

    MaterialLib::TensorParameter<PermeabilityType,
                                 IntrinsicPermeability<Matrix>, Matrix>
        K(anis_k);

    Matrix laplace =
        mock_grad_phi.trans_dphi * K.getParameterMatrix() * mock_grad_phi.dphi;

    ASSERT_NEAR(5.e-9, laplace(0, 0), 1.e-16);
    ASSERT_NEAR(5.e-9, laplace(2, 2), 1.e-16);
}
}
#endif  // OGS_USE_EIGEN
