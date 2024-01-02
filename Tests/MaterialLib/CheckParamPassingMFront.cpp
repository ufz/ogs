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

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/SolidModels/MFront/CreateMFrontGeneric.h"
#include "MaterialLib/SolidModels/MFront/Variable.h"
#include "ParameterLib/ConstantParameter.h"
#include "Tests/TestTools.h"

namespace MSM = MaterialLib::Solids::MFront;

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

    auto const strain = std::get<KV>(vars.mechanical_strain);
    auto const p = vars.liquid_phase_pressure;

    auto const stress = std::get<KV>(vars.stress);
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

auto createParams(double E, double nu, double alpha)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    auto add_param = [&parameters](const char* name, double value)
    {
        parameters.push_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>(name,
                                                                      value));
    };

    /*
     * properties: 101...199
     * inputs (gradients): X201...X299
     * outputs (forces): X301...X399
     * external state vars: X401...X499
     * X = time step (1, 2)
     */
    add_param("E", E);
    add_param("nu", nu);
    add_param("alpha", alpha);

    return parameters;
}

auto createMFront(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    auto local_coordinate_system = std::nullopt;

    const char* xml = R"XML(
        <type>MFront</type>
        <behaviour>CheckParamPassing1</behaviour>
        <library>libOgsMFrontBehaviourForUnitTests</library>
        <material_properties>
            <material_property name="YoungModulus" parameter="E"/>
            <material_property name="PoissonRatio" parameter="nu"/>
            <material_property name="ThermalExpansion" parameter="alpha"/>
        </material_properties>
        )XML";
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree config_tree(std::move(ptree), "FILENAME",
                                    &BaseLib::ConfigTree::onerror,
                                    &BaseLib::ConfigTree::onwarning);

    return MSM::createMFrontGeneric<
        3, boost::mp11::mp_list<MSM::Strain, MSM::LiquidPressure>,
        boost::mp11::mp_list<MSM::Stress, MSM::Saturation>,
        boost::mp11::mp_list<MSM::Temperature>>(
        parameters, local_coordinate_system, config_tree, false);
}

TEST(MaterialLib_CheckParamPassingMFront, SuccessTest)
{
    // Create MFront behaviour /////////////////////////////////////////////////

    auto const parameters = createParams(101, 102, 103);
    auto mfront_model = createMFront(parameters);

    // Transient data //////////////////////////////////////////////////////////

    const std::size_t num_ts = 2;
    Eigen::Vector2d ts{1, 2};
    Eigen::Vector2d Ts{1401, 2401};

    // Set initial state ///////////////////////////////////////////////////////

    auto state = mfront_model->createMaterialStateVariables();

    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    MPL::VariableArray variable_array_prev;

    // MFront thermodynamic forces
    variable_array_prev.stress = KV{1301, 1302, 1303, 1304, 1305, 1306};
    variable_array_prev.liquid_saturation = 1307;

    // MFront gradients
    variable_array_prev.mechanical_strain =
        KV{1201, 1202, 1203, 1204, 1205, 1206};
    variable_array_prev.liquid_phase_pressure = 1207;

    // MFront external state
    variable_array_prev.temperature = Ts[0];

    // Integrate material behaviour ////////////////////////////////////////////

    print_header();
    print(ts[0], variable_array_prev);

    for (std::size_t i = 1; i < num_ts; ++i)
    {
        MPL::VariableArray variable_array;

        // MFront current state
        variable_array.mechanical_strain =
            KV{2201, 2202, 2203, 2204, 2205, 2206};
        variable_array.liquid_phase_pressure = 2207;
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

        EXPECT_THAT(stress,
                    testing::Pointwise(testing::DoubleEq(),
                                       {2301, 2302, 2303, 2304, 2305, 2306}));
        EXPECT_DOUBLE_EQ(2307, S_L);
    }
}

// Consistency check that the MFront model really detects errors
TEST(MaterialLib_CheckParamPassingMFront, UnexpectedParamValues)
{
    // Create MFront behaviour /////////////////////////////////////////////////

    auto const parameters = createParams(105, 103, 102);
    auto mfront_model = createMFront(parameters);

    // Transient data //////////////////////////////////////////////////////////

    const std::size_t num_ts = 2;
    Eigen::Vector2d ts{1, 2};
    Eigen::Vector2d Ts{1401, 2401};

    // Set initial state ///////////////////////////////////////////////////////

    auto state = mfront_model->createMaterialStateVariables();

    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    MPL::VariableArray variable_array_prev;

    // MFront thermodynamic forces
    variable_array_prev.stress = KV{1301, 1302, 1303, 1304, 1305, 1306};
    variable_array_prev.liquid_saturation = 1307;

    // MFront gradients
    variable_array_prev.mechanical_strain =
        KV{1201, 1202, 1203, 1204, 1205, 1206};
    variable_array_prev.liquid_phase_pressure = 1207;

    // MFront external state
    variable_array_prev.temperature = Ts[0];

    // Integrate material behaviour ////////////////////////////////////////////

    print_header();
    print(ts[0], variable_array_prev);

    for (std::size_t i = 1; i < num_ts; ++i)
    {
        MPL::VariableArray variable_array;

        // MFront current state
        variable_array.mechanical_strain =
            KV{2201, 2202, 2203, 2204, 2205, 2206};
        variable_array.liquid_phase_pressure = 2207;
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

        EXPECT_THAT(
            stress,
            testing::Pointwise(testing::DoubleEq(),
                               // 10000 added to the first three stress
                               // components as error indicator
                               {12301, 12302, 12303, 2304, 2305, 2306}));
        EXPECT_DOUBLE_EQ(2307, S_L);
    }
}

// Consistency check that the MFront model really detects errors
TEST(MaterialLib_CheckParamPassingMFront, UnexpectedCurrentGradients)
{
    // Create MFront behaviour /////////////////////////////////////////////////

    auto const parameters = createParams(101, 102, 103);
    auto mfront_model = createMFront(parameters);

    // Transient data //////////////////////////////////////////////////////////

    const std::size_t num_ts = 2;
    Eigen::Vector2d ts{1, 2};
    Eigen::Vector2d Ts{1401, 2401};

    // Set initial state ///////////////////////////////////////////////////////

    auto state = mfront_model->createMaterialStateVariables();

    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    MPL::VariableArray variable_array_prev;

    // MFront thermodynamic forces
    variable_array_prev.stress = KV{1301, 1302, 1303, 1304, 1305, 1306};
    variable_array_prev.liquid_saturation = 1307;

    // MFront gradients
    variable_array_prev.mechanical_strain =
        KV{1201, 1202, 1203, 1204, 1205, 1206};
    variable_array_prev.liquid_phase_pressure = 1207;

    // MFront external state
    variable_array_prev.temperature = Ts[0];

    // Integrate material behaviour ////////////////////////////////////////////

    print_header();
    print(ts[0], variable_array_prev);

    for (std::size_t i = 1; i < num_ts; ++i)
    {
        MPL::VariableArray variable_array;

        // MFront current state
        variable_array.mechanical_strain =
            KV{2201 /* OK */, 2202 /* OK */, 2205 /* wrong */,
               2204 /* OK */, 2205 /* OK */, 2210 /* wrong */};
        variable_array.liquid_phase_pressure = 2208. /* wrong */;
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

        EXPECT_THAT(stress,
                    testing::Pointwise(
                        testing::DoubleEq(),
                        // 1e6 added indicating that the gradient increment is
                        // different than expected by the MFront behaviour
                        {2301, 2302, 1002303, 2304, 2305, 1002306}));
        EXPECT_DOUBLE_EQ(1002307, S_L);
    }
}

// Consistency check that the MFront model really detects errors
TEST(MaterialLib_CheckParamPassingMFront, UnexpectedPrevGradients)
{
    // Create MFront behaviour /////////////////////////////////////////////////

    auto const parameters = createParams(101, 102, 103);
    auto mfront_model = createMFront(parameters);

    // Transient data //////////////////////////////////////////////////////////

    const std::size_t num_ts = 2;
    Eigen::Vector2d ts{1, 2};
    Eigen::Vector2d Ts{1401, 2401};

    // Set initial state ///////////////////////////////////////////////////////

    auto state = mfront_model->createMaterialStateVariables();

    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    MPL::VariableArray variable_array_prev;

    // MFront thermodynamic forces
    variable_array_prev.stress = KV{1301, 1302, 1303, 1304, 1305, 1306};
    variable_array_prev.liquid_saturation = 1307;

    // MFront gradients
    variable_array_prev.mechanical_strain =
        KV{1202 /* wrong */, 1201 /* wrong */, 1203 /* OK */,
           1204 /* OK */,    1207 /* wrong */, 1206 /* OK */};
    variable_array_prev.liquid_phase_pressure = 1210;  // wrong

    // MFront external state
    variable_array_prev.temperature = Ts[0];

    // Integrate material behaviour ////////////////////////////////////////////

    print_header();
    print(ts[0], variable_array_prev);

    for (std::size_t i = 1; i < num_ts; ++i)
    {
        MPL::VariableArray variable_array;

        // MFront current state
        variable_array.mechanical_strain =
            KV{2201, 2202, 2203, 2204, 2205, 2206};
        variable_array.liquid_phase_pressure = 2207;
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

        EXPECT_THAT(stress,
                    testing::Pointwise(
                        testing::DoubleEq(),
                        // 1.1e6 added indicating that both the previous
                        // gradient and the gradient increment are different
                        // than expected by the MFront behaviour
                        {1102301, 1102302, 2303, 2304, 1102305, 2306}));
        EXPECT_DOUBLE_EQ(1102307, S_L);
    }
}

// Consistency check that the MFront model really detects errors
TEST(MaterialLib_CheckParamPassingMFront, UnexpectedPrevForces)
{
    // Create MFront behaviour /////////////////////////////////////////////////

    auto const parameters = createParams(101, 102, 103);
    auto mfront_model = createMFront(parameters);

    // Transient data //////////////////////////////////////////////////////////

    const std::size_t num_ts = 2;
    Eigen::Vector2d ts{1, 2};
    Eigen::Vector2d Ts{1401, 2401};

    // Set initial state ///////////////////////////////////////////////////////

    auto state = mfront_model->createMaterialStateVariables();

    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    MPL::VariableArray variable_array_prev;

    // MFront thermodynamic forces
    variable_array_prev.stress =
        KV{1300 /* wrong */, 1302 /* OK */,    1304 /* wrong */,
           1304 /* OK */,    1308 /* wrong */, 1306 /* OK */};
    variable_array_prev.liquid_saturation = 1300. /* wrong */;

    // MFront gradients
    variable_array_prev.mechanical_strain =
        KV{1201, 1202, 1203, 1204, 1205, 1206};
    variable_array_prev.liquid_phase_pressure = 1207;

    // MFront external state
    variable_array_prev.temperature = Ts[0];

    // Integrate material behaviour ////////////////////////////////////////////

    print_header();
    print(ts[0], variable_array_prev);

    for (std::size_t i = 1; i < num_ts; ++i)
    {
        MPL::VariableArray variable_array;

        // MFront current state
        variable_array.mechanical_strain =
            KV{2201, 2202, 2203, 2204, 2205, 2206};
        variable_array.liquid_phase_pressure = 2207;
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

        EXPECT_THAT(
            stress,
            testing::Pointwise(
                testing::DoubleEq(),
                // 1e7 added indicating that the previous value of that
                // thermodynamic force has a value not expected by the model.
                {10002301, 2302, 10002303, 2304, 10002305, 2306}));
        EXPECT_DOUBLE_EQ(10002307, S_L);
    }
}

// Consistency check that the MFront model really detects errors
TEST(MaterialLib_CheckParamPassingMFront, UnexpectedCurrentExternalState)
{
    // Create MFront behaviour /////////////////////////////////////////////////

    auto const parameters = createParams(101, 102, 103);
    auto mfront_model = createMFront(parameters);

    // Transient data //////////////////////////////////////////////////////////

    const std::size_t num_ts = 2;
    Eigen::Vector2d ts{1, 2};
    Eigen::Vector2d Ts{1401 /* OK */, 2404 /* wrong */};

    // Set initial state ///////////////////////////////////////////////////////

    auto state = mfront_model->createMaterialStateVariables();

    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    MPL::VariableArray variable_array_prev;

    // MFront thermodynamic forces
    variable_array_prev.stress = KV{1301, 1302, 1303, 1304, 1305, 1306};
    variable_array_prev.liquid_saturation = 1307;

    // MFront gradients
    variable_array_prev.mechanical_strain =
        KV{1201, 1202, 1203, 1204, 1205, 1206};
    variable_array_prev.liquid_phase_pressure = 1207;

    // MFront external state
    variable_array_prev.temperature = Ts[0];

    // Integrate material behaviour ////////////////////////////////////////////

    print_header();
    print(ts[0], variable_array_prev);

    for (std::size_t i = 1; i < num_ts; ++i)
    {
        MPL::VariableArray variable_array;

        // MFront current state
        variable_array.mechanical_strain =
            KV{2201, 2202, 2203, 2204, 2205, 2206};
        variable_array.liquid_phase_pressure = 2207;
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

        EXPECT_THAT(
            stress,
            testing::Pointwise(
                testing::DoubleEq(),
                {2301., 1e8 + 2302 /* 1e8 added because dT has an unexpected
                                     value in the MFront behaviour */
                 ,
                 2303., 2304., 2305., 2306.}));
        EXPECT_DOUBLE_EQ(2307, S_L);
    }
}

// Consistency check that the MFront model really detects errors
TEST(MaterialLib_CheckParamPassingMFront, UnexpectedPrevExternalState)
{
    // Create MFront behaviour /////////////////////////////////////////////////

    auto const parameters = createParams(101, 102, 103);
    auto mfront_model = createMFront(parameters);

    // Transient data //////////////////////////////////////////////////////////

    const std::size_t num_ts = 2;
    Eigen::Vector2d ts{1, 2};
    Eigen::Vector2d Ts{1404 /* wrong */, 2401 /* OK */};

    // Set initial state ///////////////////////////////////////////////////////

    auto state = mfront_model->createMaterialStateVariables();

    namespace MPL = MaterialPropertyLib;
    using KV = MathLib::KelvinVector::KelvinVectorType<3>;

    MPL::VariableArray variable_array_prev;

    // MFront thermodynamic forces
    variable_array_prev.stress = KV{1301, 1302, 1303, 1304, 1305, 1306};
    variable_array_prev.liquid_saturation = 1307;

    // MFront gradients
    variable_array_prev.mechanical_strain =
        KV{1201, 1202, 1203, 1204, 1205, 1206};
    variable_array_prev.liquid_phase_pressure = 1207;

    // MFront external state
    variable_array_prev.temperature = Ts[0];

    // Integrate material behaviour ////////////////////////////////////////////

    print_header();
    print(ts[0], variable_array_prev);

    for (std::size_t i = 1; i < num_ts; ++i)
    {
        MPL::VariableArray variable_array;

        // MFront current state
        variable_array.mechanical_strain =
            KV{2201, 2202, 2203, 2204, 2205, 2206};
        variable_array.liquid_phase_pressure = 2207;
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

        EXPECT_THAT(stress,
                    testing::Pointwise(
                        testing::DoubleEq(),
                        /* 1e8 added because both T and dT have unexpected
                           values in the MFront behaviour */
                        {1e8 + 2301, 1e8 + 2302, 2303., 2304., 2305., 2306.}));
        EXPECT_DOUBLE_EQ(2307, S_L);
    }
}

#endif
