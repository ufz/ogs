/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifdef OGS_USE_MFRONT

#include <gtest/gtest.h>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/SolidModels/MFront/CreateMFrontGeneric.h"
#include "MaterialLib/SolidModels/MFront/Variable.h"
#include "ParameterLib/ConstantParameter.h"
#include "Tests/TestTools.h"

static void print_header()
{
    std::cout << R"(# first column: time
# 2 column: 1th component of the gradients (StrainXX)
# 3 column: 2th component of the gradients (StrainYY)
# 4 column: 3th component of the gradients (StrainZZ)
# 5 column: 4th component of the gradients (StrainXY)
# 6 column: 5th component of the gradients (StrainXZ)
# 7 column: 6th component of the gradients (StrainYZ)
# 8 column: 7th component of the gradients (LiquidPressure)
# 9 column: 1th component of the thermodynamic forces (StressXX)
# 10 column: 2th component of the thermodynamic forces (StressYY)
# 11 column: 3th component of the thermodynamic forces (StressZZ)
# 12 column: 4th component of the thermodynamic forces (StressXY)
# 13 column: 5th component of the thermodynamic forces (StressXZ)
# 14 column: 6th component of the thermodynamic forces (StressYZ)
# 15 column: 7th component of the thermodynamic forces (Saturation)
# 16 column: stored energy
# 17 column: dissipated energy
)";
}

static void print(double const t,
                  MaterialPropertyLib::VariableArray const& vars)
{
    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    auto const strain = MaterialLib::Solids::MFront::eigenSwap45View(
                            std::get<KV>(vars.mechanical_strain))
                            .eval();
    auto const p = vars.liquid_phase_pressure;

    auto const stress =
        MaterialLib::Solids::MFront::eigenSwap45View(std::get<KV>(vars.stress))
            .eval();
    auto const S_L = vars.liquid_saturation;

    std::cout << t  //
              << ' ' << strain[0] << ' ' << strain[1] << ' ' << strain[2] << ' '
              << strain[3] << ' ' << strain[4] << ' ' << strain[5]  //
              << ' ' << p                                           //
              << ' ' << stress[0] << ' ' << stress[1] << ' ' << stress[2] << ' '
              << stress[3] << ' ' << stress[4] << ' ' << stress[5]  //
              << ' ' << S_L                                         //
              << " 0 0\n";
}

static auto createParams()
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    auto add_param = [&parameters](const char* name, double value)
    {
        parameters.push_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>(name,
                                                                      value));
    };

    add_param("E", 10e9);
    add_param("nu", 0.3);
    add_param("alpha", 1e-5);
    add_param("alpha_B", 1);
    add_param("m_chi", 1);
    add_param("S_L_res", 0);
    add_param("S_G_res", 0);
    add_param("p_b", 5e5);
    add_param("m_S", 0.4);

    return parameters;
}

static double sigma_xx_analytical(double const t)
{
    double const E = 10e9;
    double const nu = 0.3;
    double const alphaTS = 1e-5;
    double const dT = 40;
    double const p_eff = 3 * E / (3 * (1 - 2 * nu)) * alphaTS * dT;

    return -p_eff * t - 1e5;
}

static auto createMFront(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    auto local_coordinate_system = std::nullopt;

    const char* xml = R"XML(
        <type>MFront</type>
        <behaviour>ThermoPoroElasticity</behaviour>
        <material_properties>
            <material_property name="YoungModulus" parameter="E"/>
            <material_property name="PoissonRatio" parameter="nu"/>
            <material_property name="ThermalExpansion" parameter="alpha"/>
            <material_property name="BiotCoefficient" parameter="alpha_B"/>
            <material_property name="BishopsExponent" parameter="m_chi"/>
            <material_property name="ResidualLiquidSaturation" parameter="S_L_res"/>
            <material_property name="ResidualGasSaturation" parameter="S_G_res"/>
            <material_property name="BubblePressure" parameter="p_b"/>
            <material_property name="VanGenuchtenExponent_m" parameter="m_S"/>
        </material_properties>
        )XML";
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree config_tree(std::move(ptree), "FILENAME",
                                    &BaseLib::ConfigTree::onerror,
                                    &BaseLib::ConfigTree::onwarning);

    namespace MSM = MaterialLib::Solids::MFront;

    return MSM::createMFrontGeneric<
        3, boost::mp11::mp_list<MSM::Strain, MSM::LiquidPressure>,
        boost::mp11::mp_list<MSM::Stress, MSM::Saturation>,
        boost::mp11::mp_list<MSM::Temperature>>(
        parameters, local_coordinate_system, config_tree, false);
}

TEST(MaterialLib_ThermoPoroElasticityMFront, IsochoricDrainedHeating)
{
    namespace MSM = MaterialLib::Solids::MFront;

    // Create MFront behaviour /////////////////////////////////////////////////

    auto const parameters = createParams();
    auto mfront_model = createMFront(parameters);

    // Transient data //////////////////////////////////////////////////////////

    const std::size_t num_ts = 11;
    Eigen::VectorXd ts = Eigen::VectorXd::LinSpaced(num_ts, 0, 1);
    Eigen::VectorXd Ts = Eigen::VectorXd::LinSpaced(num_ts, 293, 333);

    // Set initial state ///////////////////////////////////////////////////////

    auto state = mfront_model->createMaterialStateVariables();

    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    MPL::VariableArray variable_array_prev;

    // MFront thermodynamic forces
    variable_array_prev.stress = KV{-1e5, -1e5, -1e5, 0, 0, 0};
    variable_array_prev.liquid_saturation = 1.0;

    // MFront gradients
    variable_array_prev.mechanical_strain = KV::Zero().eval();
    variable_array_prev.liquid_phase_pressure = 0.0;

    // MFront external state
    variable_array_prev.temperature = Ts[0];

    // Integrate material behaviour ////////////////////////////////////////////

    print_header();
    print(ts[0], variable_array_prev);

    for (std::size_t i = 1; i < num_ts; ++i)
    {
        MPL::VariableArray variable_array;

        // MFront current state
        variable_array.mechanical_strain = KV::Zero().eval();
        variable_array.liquid_phase_pressure = 0.0;
        variable_array.temperature = Ts[i];
        // std::cout << "T = " << Ts[i] << '\n';

        auto const t = ts[i];
        auto const dt = t - ts[i - 1];
        ParameterLib::SpatialPosition pos{};

        // Integrate
        auto solution = mfront_model->integrateStress(
            variable_array_prev, variable_array, t, pos, dt, *state);

        ASSERT_TRUE(solution);

        // Save solution for next iteration
        auto& [t_dyn_forces_data, new_state, stiffness_matrix] = *solution;

        auto const view = mfront_model->createThermodynamicForcesView();

        auto const stress = view.block(MSM::stress, t_dyn_forces_data);
        auto const S_L = view.block(MSM::saturation, t_dyn_forces_data);

        variable_array_prev = variable_array;
        variable_array_prev.stress.emplace<KV>(stress);
        variable_array_prev.liquid_saturation = S_L;

        state = std::move(new_state);

        print(t, variable_array_prev);

        EXPECT_DOUBLE_EQ(sigma_xx_analytical(t), stress[0]);
    }
}

#endif
