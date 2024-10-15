/*!
   \file
   \brief Test class for cubic law permeability model.

   \copyright
    Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 */
#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "MaterialLib/MPL/Properties/CubicLawPermeability.h"
#include "ParameterLib/ConstantParameter.h"
#include "TestMPL.h"

void testCubicLawPermeability(MPL::Property const& k_model,
                              double const aperture, bool const is_predefined_b)
{
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;
    vars.fracture_aperture = aperture;

    double const expected_k = aperture * aperture / 12.0;
    double const k = k_model.template value<double>(vars, pos, t, dt);
    ASSERT_NEAR(expected_k, k, 3.4e-15);

    double const expected_dk = is_predefined_b ? 0.0 : aperture / 6.0;
    double const dk = k_model.template dValue<double>(
        vars, MPL::Variable::fracture_aperture, pos, t, dt);
    ASSERT_NEAR(expected_dk, dk, 3.4e-15)
        << "\nExpected dk " << expected_dk << "\nCalculated dk " << dk;
}

TEST(MaterialPropertyLib, checkCubicLawPermeability)
{
    std::vector<double> const apertures = {1.0e-4, 2.e-5, 3.0e-6};
    for (double const aperture : apertures)
    {
        {  // Aperture as a variable:
            MPL::Property const& k_model = MPL::CubicLawPermeability{
                "CubicLawPermeability", nullptr /*predefined b*/};
            testCubicLawPermeability(k_model, aperture,
                                     false /*is_predefined_b*/);
        }

        {  // Aperture being pre-defined:
            ParameterLib::ConstantParameter<double> const b("b", aperture);

            MPL::Property const& k_model =
                MPL::CubicLawPermeability{"CubicLawPermeability", &b};
            testCubicLawPermeability(k_model, aperture,
                                     true /*is_predefined_b*/);
        }
    }
}
