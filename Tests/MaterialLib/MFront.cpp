/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifdef OGS_USE_MFRONT

#include "MaterialLib/SolidModels/MFront/MFront.h"

#include <gtest/gtest.h>

#include <MGIS/Behaviour/Integrate.hxx>

#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/ConstantParameter.h"

namespace MPL = MaterialPropertyLib;
using namespace MaterialLib::Solids;
template <int Dim>
using KelvinVector = MathLib::KelvinVector::KelvinVectorType<Dim>;

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
        // Parameters used by mfront model in the order of appearence in the
        // .mfront file.
        static P const young_modulus("", 1e11);
        static P const poisson_ratio("", 0.3);
        std::vector<ParameterLib::Parameter<double> const*> parameters{
            &young_modulus, &poisson_ratio};

        auto result = std::make_unique<MFront::MFront<Dim>>(
            std::move(behaviour), std::move(parameters), boost::none);
        return result;
    }
};
template <int Dim>
struct ElasticBehaviour
{
    static std::unique_ptr<MFront::MFront<Dim>> createConstitutiveRelation()
    {
        auto behaviour = mgis::behaviour::load("libOgsMFrontBehaviour",
                                               "Elasticity", hypothesis(Dim));

        using P = ParameterLib::ConstantParameter<double>;
        // Parameters used by mfront model in the order of appearence in the
        // .mfront file.
        static P const young_modulus("", 1e11);
        static P const poisson_ratio("", 0.3);
        std::vector<ParameterLib::Parameter<double> const*> parameters{
            &young_modulus, &poisson_ratio};

        auto result = std::make_unique<MFront::MFront<Dim>>(
            std::move(behaviour), std::move(parameters), boost::none);
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
        // Parameters used by mfront model in the order of appearence in the
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

        auto result = std::make_unique<MFront::MFront<Dim>>(
            std::move(behaviour), std::move(parameters), boost::none);
        return result;
    }
};

template <int Dim, typename TestBehaviour>
struct MaterialLib_SolidModelsMFront : public testing::Test
{
    MaterialLib_SolidModelsMFront()
    {
        variable_array_prev[static_cast<int>(MPL::Variable::stress)]
            .emplace<KelvinVector<Dim>>(KelvinVector<Dim>::Zero());
        variable_array_prev[static_cast<int>(MPL::Variable::strain)]
            .emplace<KelvinVector<Dim>>(KelvinVector<Dim>::Zero());
        variable_array_prev[static_cast<int>(MPL::Variable::temperature)]
            .emplace<double>(0);

        variable_array[static_cast<int>(MPL::Variable::strain)]
            .emplace<KelvinVector<Dim>>(KelvinVector<Dim>::Zero());
        variable_array[static_cast<int>(MPL::Variable::temperature)]
            .emplace<double>(0);
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
                     ElasticBehaviour<Dim>, MohrCoulombAbboSloanBehaviour<Dim>>;

TYPED_TEST_CASE(MaterialLib_SolidModelsMFront2, TestBehaviourTypes<2>);
TYPED_TEST_CASE(MaterialLib_SolidModelsMFront3, TestBehaviourTypes<3>);

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
#endif  // OGS_USE_MFRONT
