// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/SpatialPosition.h"
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

std::string thermoOsmoticProperty(std::string const& name, double const value)
{
    std::stringstream p;
    p << "    <property>\n"
         "      <name>"
      << name
      << "</name>\n"
         "      <type>Constant</type>\n"
         "      <value>"
      << value << " 0 0 " << value
      << "</value>\n"
         "    </property>\n";
    return p.str();
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
/// rejected. The name is used as the gtest parameter name.
struct RejectedMedium
{
    std::string name;
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
        evaluate({.medium_properties = thermoOsmoticProperty(
                      "thermal_osmosis_coefficient", k_T_m2_per_K_s),
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
        evaluate({.medium_properties = thermoOsmoticProperty(
                      "thermal_osmosis_permeability", epsilon_T_Pa_per_K),
                  .solid_phase_properties = ""},
                 isotropicTensor(k_intrinsic_m2), mu_liquid_Pa_s);

    expectIsotropicTensor(k_T_m2_per_K_s, result);
}

class ProcessLibThermoOsmoticCoefficientRejects
    : public ::testing::TestWithParam<RejectedMedium>
{
};

TEST_P(ProcessLibThermoOsmoticCoefficientRejects, Throws)
{
    EXPECT_ANY_THROW(evaluate(GetParam().spec,
                              isotropicTensor(k_intrinsic_m2),
                              GetParam().viscosity));
}

INSTANTIATE_TEST_SUITE_P(
    ProcessLibThermoOsmoticCoefficient,
    ProcessLibThermoOsmoticCoefficientRejects,
    ::testing::Values(
        // Defining both parametrisations at the same time is ambiguous.
        RejectedMedium{
            "BothPropertiesDefined",
            {.medium_properties =
                 thermoOsmoticProperty("thermal_osmosis_coefficient",
                                       k_T_m2_per_K_s) +
                 thermoOsmoticProperty("thermal_osmosis_permeability",
                                       epsilon_T_Pa_per_K),
             .solid_phase_properties = ""}},
        // A zero viscosity makes the conversion singular.
        RejectedMedium{"PermeabilityWithZeroViscosity",
                       {.medium_properties = thermoOsmoticProperty(
                            "thermal_osmosis_permeability", epsilon_T_Pa_per_K),
                        .solid_phase_properties = ""},
                       0.},
        // A negative viscosity flips the sign of the coefficient.
        RejectedMedium{"PermeabilityWithNegativeViscosity",
                       {.medium_properties = thermoOsmoticProperty(
                            "thermal_osmosis_permeability", epsilon_T_Pa_per_K),
                        .solid_phase_properties = ""},
                       -mu_liquid_Pa_s},
        // A leftover solid-phase thermal_osmosis_coefficient (the
        // pre-migration location of the property) must not be silently
        // ignored; it is now read from the medium, so this is rejected with a
        // migration hint instead of degrading to a silent zero.
        RejectedMedium{"LeftoverSolidPhaseCoefficient",
                       {.medium_properties = "",
                        .solid_phase_properties = thermoOsmoticProperty(
                            "thermal_osmosis_coefficient", k_T_m2_per_K_s)}}),
    [](::testing::TestParamInfo<RejectedMedium> const& info)
    { return info.param.name; });
