/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/SolidModels/MFront/MFront.h"

#include <gtest/gtest.h>

#include <MGIS/Behaviour/Integrate.hxx>

#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/ConstantParameter.h"

namespace MPL = MaterialPropertyLib;
using namespace MaterialLib::Solids;
template <int Dim>
using KelvinVector = MathLib::KelvinVector::KelvinVectorType<Dim>;
template <int Dim>
using KelvinMatrix = MathLib::KelvinVector::KelvinMatrixType<Dim>;

constexpr mgis::behaviour::Hypothesis hypothesis(int dim)
{
    switch (dim)
    {
        case 2:
            return mgis::behaviour::Hypothesis::PLANESTRAIN;
        case 3:
            return mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
    }
    // workaround for simple throw, because of gcc bugs 67371, 66026, and 80061.
    return (2 < dim && dim > 3) ? throw
                                // unreachable
                                : mgis::behaviour::Hypothesis::PLANESTRAIN;
}

template <int Dim>
struct StandardElasticityBrickBehaviour
{
    static std::unique_ptr<MFront::MFront<Dim>> createConstitutiveRelation()
    {
        auto behaviour =
            mgis::behaviour::load("libOgsMFrontBehaviour",
                                  "StandardElasticityBrick", hypothesis(Dim));

        using P = ParameterLib::ConstantParameter<double>;
        // Parameters used by mfront model in the order of appearance in the
        // .mfront file.
        static P const young_modulus("", 1e11);
        static P const poisson_ratio("", 0.3);
        std::vector<ParameterLib::Parameter<double> const*> parameters{
            &young_modulus, &poisson_ratio};
        std::map<std::string, ParameterLib::Parameter<double> const*>
            state_variables_initial_properties;

        auto result = std::make_unique<MFront::MFront<Dim>>(
            std::move(behaviour), std::move(parameters),
            std::move(state_variables_initial_properties), std::nullopt);
        return result;
    }
};

template <int Dim>
struct MohrCoulombAbboSloanBehaviour
{
    static std::unique_ptr<MFront::MFront<Dim>> createConstitutiveRelation()
    {
        auto behaviour = mgis::behaviour::load(
            "libOgsMFrontBehaviour", "MohrCoulombAbboSloan", hypothesis(Dim));

        using P = ParameterLib::ConstantParameter<double>;
        // Parameters used by mfront model in the order of appearance in the
        // .mfront file.
        static P const young_modulus("", 150e3);
        static P const poisson_ratio("", 0.3);
        static P const cohesion("", 3e1);
        static P const friction_angle("", 30);
        static P const dilatancy_angle("", 10);
        static P const transition_angle("", 29);
        static P const tension_cut_off_parameter("", 1e1);
        std::vector<ParameterLib::Parameter<double> const*> parameters{
            &young_modulus,
            &poisson_ratio,
            &cohesion,
            &friction_angle,
            &dilatancy_angle,
            &transition_angle,
            &tension_cut_off_parameter};
        std::map<std::string, ParameterLib::Parameter<double> const*>
            state_variables_initial_properties;

        auto result = std::make_unique<MFront::MFront<Dim>>(
            std::move(behaviour), std::move(parameters),
            std::move(state_variables_initial_properties), std::nullopt);
        return result;
    }
};

template <int Dim, typename TestBehaviour>
struct MaterialLib_SolidModelsMFront : public testing::Test
{
    MaterialLib_SolidModelsMFront()
    {
        variable_array_prev.stress.template emplace<KelvinVector<Dim>>(
            KelvinVector<Dim>::Zero());
        variable_array_prev.mechanical_strain
            .template emplace<KelvinVector<Dim>>(KelvinVector<Dim>::Zero());
        variable_array_prev.temperature = 0;

        variable_array.mechanical_strain.template emplace<KelvinVector<Dim>>(
            KelvinVector<Dim>::Zero());
        variable_array.temperature = 0;
        constitutive_relation = TestBehaviour::createConstitutiveRelation();
    }

    MPL::VariableArray variable_array_prev;
    MPL::VariableArray variable_array;

    double t = 0;
    ParameterLib::SpatialPosition x;
    double dt = 0;

    std::unique_ptr<MechanicsBase<Dim>> constitutive_relation;
};

template <typename TestBehaviour>
using MaterialLib_SolidModelsMFront2 =
    MaterialLib_SolidModelsMFront<2, TestBehaviour>;

template <typename TestBehaviour>
using MaterialLib_SolidModelsMFront3 =
    MaterialLib_SolidModelsMFront<3, TestBehaviour>;

template <int Dim>
using TestBehaviourTypes =
    ::testing::Types<StandardElasticityBrickBehaviour<Dim>,
                     MohrCoulombAbboSloanBehaviour<Dim>>;

TYPED_TEST_SUITE(MaterialLib_SolidModelsMFront2, TestBehaviourTypes<2>);
TYPED_TEST_SUITE(MaterialLib_SolidModelsMFront3, TestBehaviourTypes<3>);

TYPED_TEST(MaterialLib_SolidModelsMFront2, IntegrateZeroDisplacement)
{
    ASSERT_TRUE(this->constitutive_relation != nullptr);
    auto state = this->constitutive_relation->createMaterialStateVariables();

    auto solution = this->constitutive_relation->integrateStress(
        this->variable_array_prev, this->variable_array, this->t, this->x,
        this->dt, *state);

    double const epls_strain = state->getEquivalentPlasticStrain();
    double const expected_epls_strain = 0.0;
    ASSERT_LE(std::fabs(expected_epls_strain - epls_strain), 1e-10)
        << "for expected equivalent plastic strain " << expected_epls_strain
        << " and for computed equivalent plastic strain " << epls_strain;

    ASSERT_TRUE(solution != std::nullopt);
    state = std::move(std::get<1>(*solution));
    ASSERT_TRUE(state != nullptr);
    state.reset(nullptr);
    ASSERT_TRUE(state == nullptr);
}

TYPED_TEST(MaterialLib_SolidModelsMFront3, IntegrateZeroDisplacement)
{
    ASSERT_TRUE(this->constitutive_relation != nullptr);
    auto state = this->constitutive_relation->createMaterialStateVariables();

    auto solution = this->constitutive_relation->integrateStress(
        this->variable_array_prev, this->variable_array, this->t, this->x,
        this->dt, *state);

    double const epls_strain = state->getEquivalentPlasticStrain();
    double const expected_epls_strain = 0.0;
    ASSERT_LE(std::fabs(expected_epls_strain - epls_strain), 1e-10)
        << "for expected equivalent plastic strain " << expected_epls_strain
        << " and for computed equivalent plastic strain " << epls_strain;

    ASSERT_TRUE(solution != std::nullopt);
    state = std::move(std::get<1>(*solution));
    ASSERT_TRUE(state != nullptr);
    state.reset(nullptr);
    ASSERT_TRUE(state == nullptr);
}

TEST(MaterialLib_SolidModelsMFront, Conversion)
{
    using namespace MaterialLib::Solids::MFront;

    {  // vectors

        // 2D
        KelvinVector<2> const ogs2{0, 1, 2, 3};
        KelvinVector<2> const mfront2 = ogs2;

        ASSERT_EQ(mfront2, eigenSwap45View(ogs2).eval());
        ASSERT_EQ(ogs2, eigenSwap45View(mfront2).eval());

        // 3D
        KelvinVector<3> const ogs3{0, 1, 2, 3, 4, 5};
        KelvinVector<3> const mfront3{0, 1, 2, 3, 5, 4};

        ASSERT_EQ(mfront3, eigenSwap45View(ogs3).eval());
        ASSERT_EQ(ogs3, eigenSwap45View(mfront3).eval());
    }

    {  // matrices

        // 2D
        // clang-format off
        KelvinMatrix<2> ogs2;
        ogs2 <<
            0, 1, 2, 3,
            4, 5, 6, 7,
            8, 9, 10, 11,
            12, 13, 14, 15;
        KelvinMatrix<2> const mfront2 = ogs2;
        // clang-format on

        ASSERT_EQ(mfront2, eigenSwap45View(ogs2).eval());
        ASSERT_EQ(ogs2, eigenSwap45View(mfront2).eval());

        // 3D matrix
        // clang-format off
        KelvinMatrix<3> ogs3;
        ogs3 <<
            0 , 1 , 2 , 3 , 4 , 5 ,
            6 , 7 , 8 , 9 , 10, 11,
            12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29,
            30, 31, 32, 33, 34, 35;
        KelvinMatrix<3> mfront3;
        mfront3 <<                  // fourth and fifth rows/cols swapped
            0 , 1 , 2 , 3 , 5 , 4 ,
            6 , 7 , 8 , 9 , 11, 10,
            12, 13, 14, 15, 17, 16,
            18, 19, 20, 21, 23, 22,
            30, 31, 32, 33, 35, 34,
            24, 25, 26, 27, 29, 28;
        // clang-format on

        ASSERT_EQ(mfront3, eigenSwap45View(ogs3).eval());
        ASSERT_EQ(ogs3, eigenSwap45View(mfront3).eval());
    }
}
