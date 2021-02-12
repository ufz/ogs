/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <Eigen/Eigen>
#include <cmath>
#include <memory>

#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationVanGenuchten.h"
#include "MaterialLib/MPL/Properties/SaturationDependentSwelling.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ParameterLib/ConstantParameter.h"
#include "ProcessLib/RichardsMechanics/ComputeMicroPorosity.h"

using namespace ProcessLib::RichardsMechanics;
namespace MPL = MaterialPropertyLib;

namespace
{
template <int DisplacementDim>
MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>
computeElasticTangentStiffness(
    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material,
    double const t, ParameterLib::SpatialPosition const& x_position,
    double const dt, double const temperature)
{
    MPL::VariableArray variable_array;
    MPL::VariableArray variable_array_prev;

    auto const null_state = solid_material.createMaterialStateVariables();

    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    variable_array[static_cast<int>(MPL::Variable::stress)].emplace<KV>(
        KV::Zero());
    variable_array[static_cast<int>(MPL::Variable::mechanical_strain)]
        .emplace<KV>(KV::Zero());
    variable_array[static_cast<int>(MPL::Variable::temperature)]
        .emplace<double>(temperature);

    variable_array_prev[static_cast<int>(MPL::Variable::stress)].emplace<KV>(
        KV::Zero());
    variable_array_prev[static_cast<int>(MPL::Variable::mechanical_strain)]
        .emplace<KV>(KV::Zero());
    variable_array_prev[static_cast<int>(MPL::Variable::temperature)]
        .emplace<double>(temperature);

    auto&& solution = solid_material.integrateStress(
        variable_array_prev, variable_array, t, x_position, dt, *null_state);

    if (!solution)
    {
        OGS_FATAL("Computation of elastic tangent stiffness failed.");
    }

    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C =
        std::move(std::get<2>(*solution));

    return C;
}

}  // namespace

TEST(RichardsMechanics, computeMicroPorosity)
{
    static constexpr auto eps = 2e-14;
    constexpr int DisplacementDim = 2;

    static const NumLib::NewtonRaphsonSolverParameters
        nonlinear_solver_parameters{1000, 1e-8, 1e-15};

    //
    // Create properties.
    //

    // Saturation
    constexpr double residual_liquid_saturation = 0;
    constexpr double residual_gas_saturation = 0;
    constexpr double exponent = 0.4;
    constexpr double p_b = 300e6;
    MPL::Property const& saturation_micro = MPL::SaturationVanGenuchten{
        "saturation_micro", residual_liquid_saturation, residual_gas_saturation,
        exponent, p_b};

    // Swelling
    constexpr std::array swelling_pressures{1e7, 1e7, 1e7};
    constexpr std::array swelling_exponents{1.2, 1.2, 1.2};
    constexpr double lower_saturation_limit = 0;
    constexpr double upper_saturation_limit = 1;
    MPL::Property const& swelling_stress_rate =
        MPL::SaturationDependentSwelling{
            "swelling_stress_rate", swelling_pressures,     swelling_exponents,
            lower_saturation_limit, upper_saturation_limit, nullptr};

    // LinearElastic
    ParameterLib::ConstantParameter<double> E{"YoungsModulus", 150e6};
    ParameterLib::ConstantParameter<double> nu{"PoissonsRatio", 0.25};
    auto const solid_material =
        MaterialLib::Solids::LinearElasticIsotropic<DisplacementDim>{{E, nu}};

    //
    // Populate constants
    //
    ParameterLib::SpatialPosition const pos;
    constexpr double t0 = 0;
    constexpr double t_max = 200e3;
    double const dt = 100;
    double const T = 293.15;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>::
        identity2;
    auto const C_el =
        computeElasticTangentStiffness(solid_material, t0, pos, dt, T);
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const&
        I_2_C_el_inverse = identity2.transpose() * C_el.inverse();
    double const rho_LR_m = 1e3;
    double const mu_LR = 1e-3;
    double const alpha_bar = 1e-14;
    double const alpha_B = 1;
    double const phi_M = 0.45;

    auto saturation = [&](double const p_L) {
        MPL::VariableArray v;
        v[static_cast<int>(MPL::Variable::capillary_pressure)] = -p_L;
        return saturation_micro.template value<double>(v, pos, t0, dt);
    };

    //
    // Drive the computation by macro pressure. Saturation, steady, and
    // desaturation.
    //

    // The time coordinates are simultaneously also the test points.
    std::vector const t_coords{
        // saturation is quick to rich maximum,
        0 * t_max, 0.05 * t_max, 0.2 * t_max, 0.3 * t_max,
        // desaturation takes much longer.
        0.35 * t_max, 0.5 * t_max, t_max};
    std::vector const pressure_M_values{// saturation
                                        -200e6, -200e6, 0., 0.,
                                        // desaturation
                                        0., -200e6, -200e6};
    MathLib::PiecewiseLinearInterpolation const pressure_M{t_coords,
                                                           pressure_M_values};

    //
    // Initial state
    //

    MicroPorosityStateSpace<DisplacementDim> state_prev = {
        0.45 - 0.25,  // phi_m
        0,            // e_sw
        -200e6,       // p_L_m
        {0, 0, 0, 0}  // sigma_sw
    };

    auto state = state_prev;

    std::vector<MicroPorosityStateSpace<DisplacementDim>> results;
    for (double t = t0 + dt; t <= t_max; t += dt)
    {
        double const p_L = pressure_M.getValue(t);
        double const S_L_m_prev = saturation(state_prev.p_L_m);
        //auto const [delta_phi_m, delta_e_sw, delta_sigma_sw, delta_p_L_m] =
        auto const state_increment = computeMicroPorosity<DisplacementDim>(
            I_2_C_el_inverse,
            rho_LR_m,  // for simplification equal to rho_LR_M
            mu_LR, alpha_bar, alpha_B, phi_M, p_L, state_prev.p_L_m,
            MaterialPropertyLib::VariableArray{}, S_L_m_prev, state_prev.phi_m,
            pos, t, dt, saturation_micro, swelling_stress_rate,
            nonlinear_solver_parameters);

        // push back state
        state_prev = state;

        // update
        state += state_increment;

        if (std::find_if(begin(t_coords), end(t_coords), [&](auto const value) {
                return std::abs(t - value) < eps;
            }) != end(t_coords))
        {
            results.push_back(state);
            /* Keep for possible result updates
            std::cout << std::setprecision(17)
                      << "MicroPorosityStateSpace<DisplacementDim>{ "
                      << state.phi_m << ", " << state.e_sw << ", "
                      << state.p_L_m << ", {" << state.sigma_sw[0] << ", "
                      << state.sigma_sw[1] << ", " << state.sigma_sw[2] << ", "
                      << state.sigma_sw[3] << " }}\n";
            */
        }
    }

    //
    // Expected results
    //
    std::vector expected_results = {
        MicroPorosityStateSpace<DisplacementDim>{
            0.20000000000000001, 0, -200000000, {0, 0, 0, 0}},
        MicroPorosityStateSpace<DisplacementDim>{
            0.20678794335787565,
            0.012341715196137344,
            -87023875.75987792,
            {-1234171.5196137347, -1234171.5196137347, -1234171.5196137347, 0}},
        MicroPorosityStateSpace<DisplacementDim>{
            0.20984996711734061,
            0.017909031122437243,
            -4693546.2948520193,
            {-1790903.1122437208, -1790903.1122437208, -1790903.1122437208, 0}},
        MicroPorosityStateSpace<DisplacementDim>{
            0.20987701765630837,
            0.017958213920560351,
            23234.984113500199,
            {-1795821.3920560319, -1795821.3920560319, -1795821.3920560319, 0}},
        MicroPorosityStateSpace<DisplacementDim>{
            0.20475472126658667,
            0.0086449477574298151,
            -123250534.84079438,
            {-864494.77574297995, -864494.77574297995, -864494.77574297995, 0}},
        MicroPorosityStateSpace<DisplacementDim>{
            0.20005534659406693,
            0.0001006301710302317,
            -199863745.41244847,
            {-10063.017103021764, -10063.017103021764, -10063.017103021764,
             0}}};

    auto eps_equal = [](double const a, double const b) {
        return std::abs(a - b) < eps;
    };

    auto const [mismatch_it_expected, mismatch_it_results] = std::mismatch(
        begin(expected_results), end(expected_results), begin(results),
        [&](auto const& a, auto const& b) {
            EXPECT_TRUE(eps_equal(a.phi_m, b.phi_m)) << "with eps = " << eps;
            EXPECT_TRUE(eps_equal(a.e_sw, b.e_sw)) << "with eps = " << eps;
            EXPECT_TRUE(eps_equal(a.p_L_m / 1e9, b.p_L_m / 1e9))
                << "with eps = " << eps;
            EXPECT_TRUE((a.sigma_sw - b.sigma_sw).norm() / 1e9 < eps)
                << "with eps = " << eps;
            // Divide pressures and stresses by GPa to get approxmately
            // O(1) numbers.
            return eps_equal(a.phi_m, b.phi_m) && eps_equal(a.e_sw, b.e_sw) &&
                   eps_equal(a.p_L_m / 1e9, b.p_L_m / 1e9) &&
                   (a.sigma_sw - b.sigma_sw).norm() / 1e9 < eps;
        });

    ASSERT_TRUE(mismatch_it_expected == end(expected_results))
        << "Mismatch for expected result\n\t" << *mismatch_it_expected
        << "\n\tfound. Corresponding result is\n\t" << *mismatch_it_results;
    ASSERT_TRUE(mismatch_it_results == end(results));
}
