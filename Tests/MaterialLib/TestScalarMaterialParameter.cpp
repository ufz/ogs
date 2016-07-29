/*!
   \file  TestScalarMaterialParameter.cpp
   \brief Test classes for scalar material parameters.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <algorithm>    // std::max

#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/ScalarParameter.h"
#include "MaterialLib/ConstantScalarModel.h"

#include "MaterialLib/Solid/Density/LinearSolidDensityModel.h"
#include "MaterialLib/Fluid/Density/IdealGasLaw.h"

#include "MaterialLib/Fluid/Viscosity/PressureDependentViscosity.h"

namespace
{

// For local assembly
template<typename T_MAT>class TestScalar
{
    public:
    TestScalar(const T_MAT &mat) : _mat(const_cast<T_MAT*>(&mat))
    {
    }

    // Test function
    template<typename... Args> double getMatParameterTest(Args... args) const
    {
          return _mat->getValue(args...);
    }
    private:
      T_MAT *_mat;
};

using namespace MaterialLib;

TEST(Material, checkDensity)
{
    //-- Solid --------------------------------------------------------------
    // Constant
    const double rho = 2080.;
    MaterialLib::ScalarParameter<DensityType, ConstantScalarModel> s_density(rho);
    ASSERT_NEAR(rho, s_density.getValue(), 1.e-10);
    ASSERT_EQ(DensityType::CONSTANT, s_density.getType());

    // Linear
    MaterialLib::ScalarParameter<DensityType, LinearSolidDensityModel>
    lin_density(20., 1000., 100., 900.);
    const double lin_den_expected = (900. - 1000.)/(100. - 20.)
                                    * (50.-20.) + 1000.;
    ASSERT_NEAR(lin_den_expected, lin_density.getValue(50.), 1.e-10);

    //-- Fluid --------------------------------------------------------------
    const double molar_air = 28.96;
    MaterialLib::ScalarParameter<DensityType, IdealGasLaw> air_density(molar_air);
    const double T = 290.;
    const double p = 1.e+5;
    const double R = PhysicalConstant::IdealGasConstant;
    const double expected_air_dens = molar_air * p /(R * T);
    ASSERT_NEAR(expected_air_dens, air_density.getValue(T, p), 1.e-10);

    const double expected_d_air_dens_dT = -molar_air * p /(R * T * T);
    ASSERT_NEAR(expected_d_air_dens_dT, air_density.getdValue(T, p, 0), 1.e-10);

    const double expected_d_air_dens_dp = molar_air /(R * T);
    ASSERT_NEAR(expected_d_air_dens_dp, air_density.getdValue(T, p, 1), 1.e-10);

    TestScalar< MaterialLib::ScalarParameter<DensityType, IdealGasLaw> > test0(air_density);
    ASSERT_NEAR(expected_air_dens, test0.getMatParameterTest(T, p), 1.e-10);

    // Check polymorphy
    MaterialLib::ParameterBase<DensityType> *den_ptr_base = &air_density;
    if(den_ptr_base->getType() == DensityType::IDEAL_GAS)
    {
        MaterialLib::ScalarParameter<DensityType, IdealGasLaw> *den_ptr
            = static_cast< MaterialLib::ScalarParameter<DensityType, IdealGasLaw>* >(den_ptr_base);
        ASSERT_NEAR(molar_air * 2.e+5 /(R * 253.), den_ptr->getValue(253., 2.e+5), 1.e-10);
    }
}

TEST(Material, checkViscosity)
{
    // Constant
    const double mu = 1.e-3;
    MaterialLib::ScalarParameter<ViscosityType, ConstantScalarModel> mu_const(mu);
    ASSERT_NEAR(mu, mu_const.getValue(), 1.e-10);
    ASSERT_EQ(ViscosityType::CONSTANT, mu_const.getType());

    // Pressure dependent
    const double mu0 = 1.e-3;
    const double p0 = 1.e+5;
    const double gamma = 10.;
    MaterialLib::ScalarParameter<ViscosityType, PressureDependentViscosity>
        mu_pre(mu0, p0, gamma);
    const double mu_pre_expected = mu0 * (1 + gamma* (std::max(1.e+6, 0.) - p0) );
    ASSERT_NEAR(mu_pre_expected, mu_pre.getValue(1.e+6), 1.e-10);
}

}
