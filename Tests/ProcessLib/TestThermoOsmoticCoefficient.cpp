// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/Common/ThermoOsmosis/CheckThermoOsmosisProperties.h"
#include "ProcessLib/Common/ThermoOsmosis/ThermoOsmoticCoefficient.h"
#include "Tests/MaterialLib/TestMPL.h"

namespace
{
// Values taken from the ThermoOsmosis/Column benchmarks, where the
// thermal_osmosis_coefficient and thermal_osmosis_permeability
// parametrisations are equivalent:
// epsilon_T k / mu = 5400 * 5e-17 / 1e-3 = 2.7e-10.
constexpr double k_intrinsic_m2 = 5e-17;
constexpr double mu_liquid_Pa_s = 1e-3;
constexpr double epsilon_T_Pa_per_K = 5400;
constexpr double k_T_m2_per_K_s = 2.7e-10;

/// Medium to assemble: the given medium-level and solid-phase properties. The
/// intrinsic permeability and the liquid viscosity are not read from the
/// medium; they are passed to the helper by its caller.
struct MediumSpec
{
    std::string medium_properties;
    std::string solid_phase_properties;
};

/// A Constant property with the given components, which are read as a scalar
/// for a single component and as a row-major tensor otherwise, following
/// MaterialPropertyLib::fromVector().
std::string constantProperty(std::string const& name,
                             std::vector<double> const& components)
{
    std::stringstream p;
    p << "    <property>\n"
         "      <name>"
      << name
      << "</name>\n"
         "      <type>Constant</type>\n"
         "      <value>";
    for (auto const component : components)
    {
        p << component << " ";
    }
    p << "</value>\n"
         "    </property>\n";
    return p.str();
}

/// thermal_osmosis_permeability is a scalar property.
std::string permeabilityProperty(double const epsilon_T)
{
    return constantProperty("thermal_osmosis_permeability", {epsilon_T});
}

/// thermal_osmosis_coefficient is a tensor property; the isotropic tensor with
/// the given diagonal entry.
std::string coefficientProperty(double const k_T)
{
    return constantProperty("thermal_osmosis_coefficient", {k_T, 0., 0., k_T});
}

std::string makeMedium(MediumSpec const& spec)
{
    std::stringstream m;
    m << "<medium>\n";
    if (!spec.solid_phase_properties.empty())
    {
        m << "  <phases>\n"
             "    <phase>\n"
             "      <type>Solid</type>\n"
             "      <properties>\n"
          << spec.solid_phase_properties
          << "      </properties>\n"
             "    </phase>\n"
             "  </phases>\n";
    }
    // The porosity is never read by the helper. It keeps <properties>
    // non-empty for the case without any thermo-osmosis property, which
    // createMedium() rejects as a medium with neither phases nor properties.
    m << "  <properties>\n"
         "    <property>\n"
         "      <name>porosity</name>\n"
         "      <type>Constant</type>\n"
         "      <value>0.2</value>\n"
         "    </property>\n"
      << spec.medium_properties
      << "  </properties>\n"
         "</medium>\n";
    return m.str();
}

Eigen::Matrix<double, 2, 2> isotropicTensor(double const diagonal_entry)
{
    Eigen::Matrix<double, 2, 2> tensor = Eigen::Matrix<double, 2, 2>::Zero();
    tensor(0, 0) = diagonal_entry;
    tensor(1, 1) = diagonal_entry;
    return tensor;
}

Eigen::Matrix<double, 2, 2> evaluate(
    MediumSpec const& spec,
    Eigen::Matrix<double, 2, 2> const& intrinsic_permeability,
    double const liquid_dynamic_viscosity)
{
    auto const medium = Tests::createTestMaterial(makeMedium(spec), 2);
    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    return ProcessLib::getThermoOsmoticCoefficient<2>(
        *medium, variable_array, pos,
        /*t=*/0., /*dt=*/1., intrinsic_permeability, liquid_dynamic_viscosity);
}

/// Checks that the tensor is the isotropic one with the given diagonal, with
/// exactly zero off-diagonal entries.
void expectIsotropicTensor(double const expected_diagonal_entry,
                           Eigen::Matrix<double, 2, 2> const& result)
{
    EXPECT_NEAR(expected_diagonal_entry, result(0, 0), 1e-24);
    EXPECT_NEAR(expected_diagonal_entry, result(1, 1), 1e-24);
    EXPECT_EQ(0., result(0, 1));
    EXPECT_EQ(0., result(1, 0));
}

/// A medium the coefficient cannot be evaluated for, and the reason it is
/// rejected. The name is used as the gtest parameter name;
/// expected_message_fragment identifies which of the helper's fatal errors is
/// expected, so that a case cannot pass on an unrelated failure.
struct RejectedMedium
{
    std::string name;
    std::string expected_message_fragment;
    MediumSpec spec;
    double viscosity = mu_liquid_Pa_s;
};

std::ostream& operator<<(std::ostream& os, RejectedMedium const& c)
{
    return os << c.name;
}
}  // namespace

// Without either property there is no thermo-osmosis. The helper returns an
// exact Eigen::Matrix::Zero(), so the test pins exact zeros rather than
// isZero()'s default tolerance.
TEST(ProcessLibThermoOsmoticCoefficient, NoPropertyYieldsZeroTensor)
{
    auto const result =
        evaluate({.medium_properties = "", .solid_phase_properties = ""},
                 isotropicTensor(k_intrinsic_m2), mu_liquid_Pa_s);

    EXPECT_TRUE((result.array() == 0.).all()) << "got\n" << result;
}

// thermal_osmosis_coefficient is used as given, and neither the permeability
// nor the viscosity the caller passes in enters the result.
TEST(ProcessLibThermoOsmoticCoefficient, CoefficientIsPassedThrough)
{
    auto const result =
        evaluate({.medium_properties = coefficientProperty(k_T_m2_per_K_s),
                  .solid_phase_properties = ""},
                 isotropicTensor(2. * k_intrinsic_m2), 2. * mu_liquid_Pa_s);

    expectIsotropicTensor(k_T_m2_per_K_s, result);
}

// thermal_osmosis_permeability is converted via k_T = epsilon_T k / mu, and the
// conversion reproduces the equivalent thermal_osmosis_coefficient to
// floating-point precision (1e-24 absolute, i.e. ~1e-14 relative).
// This pins the two parametrisations against each other at the same k and mu
// the caller uses for the Darcy term.
TEST(ProcessLibThermoOsmoticCoefficient, PermeabilityIsConvertedToCoefficient)
{
    auto const result =
        evaluate({.medium_properties = permeabilityProperty(epsilon_T_Pa_per_K),
                  .solid_phase_properties = ""},
                 isotropicTensor(k_intrinsic_m2), mu_liquid_Pa_s);

    expectIsotropicTensor(k_T_m2_per_K_s, result);
}

// epsilon_T is a scalar, so k_T is k scaled by epsilon_T / mu and inherits k's
// anisotropy, including its off-diagonal entries. With
// epsilon_T / mu = 5400 / 1e-3 = 5.4e6 the expected entries are k's, scaled:
//   k_T = 5.4e6 * [[5e-17, 0], [2e-17, 5e-17]]
//       = [[2.7e-10, 0], [1.08e-10, 2.7e-10]]
TEST(ProcessLibThermoOsmoticCoefficient, ScalarPermeabilityScalesAnisotropicK)
{
    Eigen::Matrix<double, 2, 2> k;
    k << k_intrinsic_m2, 0., 2e-17, k_intrinsic_m2;

    auto const result =
        evaluate({.medium_properties = permeabilityProperty(epsilon_T_Pa_per_K),
                  .solid_phase_properties = ""},
                 k, mu_liquid_Pa_s);

    EXPECT_NEAR(2.7e-10, result(0, 0), 1e-24);
    EXPECT_EQ(0., result(0, 1));
    EXPECT_NEAR(1.08e-10, result(1, 0), 1e-24);
    EXPECT_NEAR(2.7e-10, result(1, 1), 1e-24);
}

// A tensor-valued thermal_osmosis_permeability is not a valid parametrisation:
// epsilon_T is a scalar, and the value access rejects the tensor rather than
// taking a component of it.
TEST(ProcessLibThermoOsmoticCoefficient, TensorPermeabilityIsRejected)
{
    MediumSpec const spec{.medium_properties = constantProperty(
                              "thermal_osmosis_permeability",
                              {epsilon_T_Pa_per_K, 0., 0., epsilon_T_Pa_per_K}),
                          .solid_phase_properties = ""};

    try
    {
        evaluate(spec, isotropicTensor(k_intrinsic_m2), mu_liquid_Pa_s);
        FAIL() << "expected a fatal error about the requested type";
    }
    catch (std::runtime_error const& e)
    {
        EXPECT_THAT(
            e.what(),
            ::testing::HasSubstr(
                "'thermal_osmosis_permeability' defined for medium 0 is not "
                "of the requested type 'double' but a 2x2-matrix"));
    }
}

class ProcessLibThermoOsmoticCoefficientRejects
    : public ::testing::TestWithParam<RejectedMedium>
{
};

TEST_P(ProcessLibThermoOsmoticCoefficientRejects, Throws)
{
    auto const& rejected_medium = GetParam();
    try
    {
        evaluate(rejected_medium.spec, isotropicTensor(k_intrinsic_m2),
                 rejected_medium.viscosity);
        FAIL() << "expected a fatal error mentioning \""
               << rejected_medium.expected_message_fragment << '"';
    }
    catch (std::runtime_error const& e)
    {
        EXPECT_THAT(
            e.what(),
            ::testing::HasSubstr(rejected_medium.expected_message_fragment));
    }
}

INSTANTIATE_TEST_SUITE_P(
    ProcessLibThermoOsmoticCoefficient,
    ProcessLibThermoOsmoticCoefficientRejects,
    ::testing::Values(
        // A zero viscosity makes the conversion singular.
        RejectedMedium{
            "PermeabilityWithZeroViscosity",
            "viscosity must be > 0",
            {.medium_properties = permeabilityProperty(epsilon_T_Pa_per_K),
             .solid_phase_properties = ""},
            0.},
        // A negative viscosity flips the sign of the coefficient.
        RejectedMedium{
            "PermeabilityWithNegativeViscosity",
            "viscosity must be > 0",
            {.medium_properties = permeabilityProperty(epsilon_T_Pa_per_K),
             .solid_phase_properties = ""},
            -mu_liquid_Pa_s}),
    [](::testing::TestParamInfo<RejectedMedium> const& info)
    { return info.param.name; });

// The parametrisation itself is checked once at process creation, by
// checkThermoOsmosisProperties(), not on every coefficient evaluation.
class ProcessLibCheckThermoOsmosisPropertiesRejects
    : public ::testing::TestWithParam<RejectedMedium>
{
};

TEST_P(ProcessLibCheckThermoOsmosisPropertiesRejects, Throws)
{
    auto const& rejected_medium = GetParam();
    auto const medium =
        Tests::createTestMaterial(makeMedium(rejected_medium.spec), 2);
    try
    {
        ProcessLib::checkThermoOsmosisProperties(*medium);
        FAIL() << "expected a fatal error mentioning \""
               << rejected_medium.expected_message_fragment << '"';
    }
    catch (std::runtime_error const& e)
    {
        EXPECT_THAT(
            e.what(),
            ::testing::HasSubstr(rejected_medium.expected_message_fragment));
    }
}

INSTANTIATE_TEST_SUITE_P(
    ProcessLibCheckThermoOsmosisProperties,
    ProcessLibCheckThermoOsmosisPropertiesRejects,
    ::testing::Values(
        // Defining both parametrisations at the same time is ambiguous.
        RejectedMedium{
            "BothPropertiesDefined",
            "cannot be defined at the same time",
            {.medium_properties = coefficientProperty(k_T_m2_per_K_s) +
                                  permeabilityProperty(epsilon_T_Pa_per_K),
             .solid_phase_properties = ""}},
        // A leftover solid-phase property (thermal_osmosis_coefficient's
        // pre-migration location) must not be silently ignored; both
        // properties are read from the medium, so either one on the solid
        // phase is rejected with a migration hint instead of degrading to a
        // silent zero.
        RejectedMedium{
            "LeftoverSolidPhaseCoefficient",
            "thermal_osmosis_coefficient is defined on the solid "
            "phase",
            {.medium_properties = "",
             .solid_phase_properties = coefficientProperty(k_T_m2_per_K_s)}},
        RejectedMedium{"SolidPhasePermeability",
                       "thermal_osmosis_permeability is defined on the solid "
                       "phase",
                       {.medium_properties = "",
                        .solid_phase_properties =
                            permeabilityProperty(epsilon_T_Pa_per_K)}}),
    [](::testing::TestParamInfo<RejectedMedium> const& info)
    { return info.param.name; });

// A medium the coefficient can be evaluated for passes the check. Both
// parametrisations, and a medium without any thermo-osmosis property, are
// accepted.
TEST(ProcessLibCheckThermoOsmosisProperties, AcceptsValidParametrisations)
{
    for (auto const& spec :
         {MediumSpec{.medium_properties = "", .solid_phase_properties = ""},
          MediumSpec{.medium_properties = coefficientProperty(k_T_m2_per_K_s),
                     .solid_phase_properties = ""},
          MediumSpec{
              .medium_properties = permeabilityProperty(epsilon_T_Pa_per_K),
              .solid_phase_properties = ""}})
    {
        auto const medium = Tests::createTestMaterial(makeMedium(spec), 2);
        EXPECT_NO_THROW(ProcessLib::checkThermoOsmosisProperties(*medium));
    }
}
