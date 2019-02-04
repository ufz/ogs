/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifdef OGS_USE_MFRONT

#include <gtest/gtest.h>
#include <MGIS/Behaviour/Integrate.hxx>

#include "MaterialLib/SolidModels/MFront/MFront.h"
#include "ParameterLib/ConstantParameter.h"

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
struct ElasticBehaviour
{
    static std::unique_ptr<MFront::MFront<Dim>> createConstitutiveRelation()
    {
        auto behaviour = mgis::behaviour::load("libOgsMFrontBehaviour.so",
                                               "Elasticity", hypothesis(Dim));

        using P = ProcessLib::ConstantParameter<double>;
        // Parameters used by mfront model in the order of appearence in the
        // .mfront file.
        static P const young_modulus("", 1e11);
        static P const poisson_ratio("", 1e9);
        std::vector<ProcessLib::Parameter<double> const*> parameters{
            &young_modulus, &poisson_ratio};

        auto result = std::make_unique<MFront::MFront<Dim>>(
            std::move(behaviour), std::move(parameters));
        return result;
    }
};

template <int Dim, typename TestBehaviour>
struct MaterialLib_SolidModelsMFront : public testing::Test
{
    MaterialLib_SolidModelsMFront()
    {
        constitutive_relation = TestBehaviour::createConstitutiveRelation();
    }

    KelvinVector<Dim> const eps_prev = KelvinVector<Dim>::Zero();
    KelvinVector<Dim> const eps = KelvinVector<Dim>::Zero();
    KelvinVector<Dim> const sigma_prev = KelvinVector<Dim>::Zero();

    double t = 0;
    ProcessLib::SpatialPosition x;
    double dt = 0;
    double T = 0;

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
    ::testing::Types<ElasticBehaviour<Dim>>;

TYPED_TEST_CASE(MaterialLib_SolidModelsMFront2, TestBehaviourTypes<2>);
TYPED_TEST_CASE(MaterialLib_SolidModelsMFront3, TestBehaviourTypes<3>);

TYPED_TEST(MaterialLib_SolidModelsMFront2, IntegrateZeroDisplacement)
{
    ASSERT_TRUE(this->constitutive_relation != nullptr);
    auto state = this->constitutive_relation->createMaterialStateVariables();

    auto solution = this->constitutive_relation->integrateStress(
        this->t, this->x, this->dt, this->eps_prev, this->eps, this->sigma_prev,
        *state, this->T);

    ASSERT_TRUE(solution != boost::none);
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
        this->t, this->x, this->dt, this->eps_prev, this->eps, this->sigma_prev,
        *state, this->T);

    ASSERT_TRUE(solution != boost::none);
    state = std::move(std::get<1>(*solution));
    ASSERT_TRUE(state != nullptr);
    state.reset(nullptr);
    ASSERT_TRUE(state == nullptr);
}
#endif  // OGS_USE_MFRONT
