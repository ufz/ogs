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
#include <spdlog/fmt/bundled/format.h>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/SolidModels/MFront/CreateMFrontGeneric.h"
#include "MaterialLib/SolidModels/MFront/Variable.h"
#include "Tests/TestTools.h"

namespace MSM = MaterialLib::Solids::MFront;

template <typename ExtStateVars>
auto createMFront(std::string const& behaviour)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    auto local_coordinate_system = std::nullopt;

    auto const xml = fmt::format(R"XML(
        <type>MFront</type>
        <behaviour>{}</behaviour>
        <library>libOgsMFrontBehaviourForUnitTests</library>
        <material_properties />
        )XML",
                                 behaviour);
    auto ptree = Tests::readXml(xml.c_str());
    BaseLib::ConfigTree config_tree(std::move(ptree), "FILENAME",
                                    &BaseLib::ConfigTree::onerror,
                                    &BaseLib::ConfigTree::onwarning);

    return MSM::createMFrontGeneric<
        3, boost::mp11::mp_list<MSM::Strain, MSM::LiquidPressure>,
        boost::mp11::mp_list<MSM::Stress, MSM::Saturation>, ExtStateVars>(
        parameters, local_coordinate_system, config_tree, false);
}

struct TestDataBase
{
    virtual Eigen::MatrixXd dsig_deps() const
    {
        Eigen::MatrixXd dsig_deps(6, 6);

        for (Eigen::Index r = 0; r < 6; ++r)
        {
            for (Eigen::Index c = 0; c < 6; ++c)
            {
                dsig_deps(r, c) = 1000 + 6 * r + c + 1;
            }
        }

        return dsig_deps;
    }

    virtual Eigen::VectorXd dsig_dp() const
    {
        Eigen::VectorXd dsig_dp(6);
        dsig_dp << 2001, 2002, 2003, 2004, 2005, 2006;
        return dsig_dp;
    }

    virtual Eigen::MatrixXd dsat_deps() const
    {
        Eigen::MatrixXd dsat_deps(1, 6);
        dsat_deps << 3001, 3002, 3003, 3004, 3005, 3006;
        return dsat_deps;
    }

    virtual double dsat_dp() const { return 4001; }
};

class MaterialLib_CheckStiffnessMatrixMFront : public ::testing::Test
{
protected:
    template <typename ExtStateVars = boost::mp11::mp_list<MSM::Temperature>>
    void run(std::string const& behaviour, TestDataBase const& data) const
    {
        run<ExtStateVars>(behaviour, data,
                          [](auto const& /*view*/, auto const& /*data*/) {});
    }

    template <typename ExtStateVars, typename AdditionalChecks>
    void run(std::string const& behaviour, TestDataBase const& data,
             AdditionalChecks&& checks) const
    {
        // Create MFront behaviour
        // /////////////////////////////////////////////////

        auto mfront_model = createMFront<ExtStateVars>(behaviour);

        // Transient data
        // //////////////////////////////////////////////////////////

        const std::size_t num_ts = 2;
        Eigen::Vector2d ts{1, 2};
        Eigen::Vector2d Ts{1401, 2401};

        // Set initial state
        // ///////////////////////////////////////////////////////

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

        // Integrate material behaviour
        // ////////////////////////////////////////////

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
            auto& [t_dyn_forces_data, new_state, tangent_operator_data] =
                *solution;

            auto const view = mfront_model->createThermodynamicForcesView();

            auto const stress = view.block(MSM::stress, t_dyn_forces_data);
            auto const S_L = view.block(MSM::saturation, t_dyn_forces_data);

            variable_array_prev = variable_array;
            variable_array_prev.stress.emplace<KV>(stress);
            variable_array_prev.liquid_saturation = S_L;

            state = std::move(new_state);

            EXPECT_THAT(
                stress,
                testing::Pointwise(testing::DoubleEq(),
                                   {2301, 2302, 2303, 2304, 2305, 2306}));
            EXPECT_DOUBLE_EQ(2307, S_L);

            auto const blocks_view =
                mfront_model->createTangentOperatorBlocksView();

            EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, data.dsig_deps(),
                                (blocks_view.block(MSM::stress, MSM::strain,
                                                   tangent_operator_data)));

            EXPECT_PRED_FORMAT2(
                Tests::EigenIsNear{}, data.dsig_dp(),
                (blocks_view.block(MSM::stress, MSM::liquid_pressure,
                                   tangent_operator_data)));

            EXPECT_PRED_FORMAT2(Tests::EigenIsNear{}, data.dsat_deps(),
                                (blocks_view.block(MSM::saturation, MSM::strain,
                                                   tangent_operator_data)));

            EXPECT_DOUBLE_EQ(
                data.dsat_dp(),
                (blocks_view.block(MSM::saturation, MSM::liquid_pressure,
                                   tangent_operator_data)));

            checks(blocks_view, tangent_operator_data);
        }
    }
};

class AllBlocksProvidedTest : public MaterialLib_CheckStiffnessMatrixMFront,
                              public testing::WithParamInterface<const char*>
{
};

INSTANTIATE_TEST_SUITE_P(
    MaterialLib_CheckStiffnessMatrixMFront_Parameterized,
    AllBlocksProvidedTest,
    testing::Values(
        "CheckStiffnessMatrixAllBlocks",
        // Reordering should make no difference
        "CheckStiffnessMatrixAllBlocksReordered",
        // Excess blocks are ignored and should also make no difference
        "CheckStiffnessMatrixExcessBlocks"));

TEST_P(AllBlocksProvidedTest, Test)
{
    run(GetParam(), TestDataBase{});
}

struct TestDataMissingBlocks : TestDataBase
{
    Eigen::MatrixXd dsig_deps() const override
    {
        return Eigen::MatrixXd::Zero(6, 6);
    }

    double dsat_dp() const override { return 0; }
};

TEST_F(MaterialLib_CheckStiffnessMatrixMFront, SomeBlocksMissing)
{
    run("CheckStiffnessMatrixMissingBlocks", TestDataMissingBlocks{});
}

TEST_F(MaterialLib_CheckStiffnessMatrixMFront, DSigmaDTBlock)
{
    run<boost::mp11::mp_list<MSM::Temperature>>(
        "CheckStiffnessMatrixExcessBlocks", TestDataBase{},
        [](auto const& blocks_view, auto const& tangent_operator_data)
        {
            Eigen::VectorXd dsig_dT(6);
            dsig_dT << 5001, 5002, 5003, 5004, 5005, 5006;

            EXPECT_PRED_FORMAT2(
                Tests::EigenIsNear{}, dsig_dT,
                (blocks_view.block(MSM::stress, MSM::temperature,
                                   tangent_operator_data)));
        });
}

#endif
