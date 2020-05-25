/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cassert>


#include "MathLib/Nonlinear/Root1D.h"
#include "MathLib/ODE/ODESolverBuilder.h"

#include "MaterialLib/Adsorption/Adsorption.h"
#include "MaterialLib/Adsorption/ReactionInert.h"
#include "MaterialLib/Adsorption/ReactionSinusoidal.h"

#include "TESLocalAssemblerInner.h"
#include "TESReactionAdaptor.h"

namespace ProcessLib
{
namespace TES
{
std::unique_ptr<TESFEMReactionAdaptor> TESFEMReactionAdaptor::newInstance(
    TESLocalAssemblerData const& data)
{
    auto const* ads = data.ap.react_sys.get();
    if (dynamic_cast<Adsorption::AdsorptionReaction const*>(ads) != nullptr)
    {
        return std::make_unique<TESFEMReactionAdaptorAdsorption>(data);
    }
    if (dynamic_cast<Adsorption::ReactionInert const*>(ads) != nullptr)
    {
        return std::make_unique<TESFEMReactionAdaptorInert>(data);
    }
    if (dynamic_cast<Adsorption::ReactionSinusoidal const*>(ads) != nullptr)
    {
        return std::make_unique<TESFEMReactionAdaptorSinusoidal>(data);
    }
    if (dynamic_cast<Adsorption::ReactionCaOH2 const*>(ads) != nullptr)
    {
        return std::make_unique<TESFEMReactionAdaptorCaOH2>(data);
    }

    OGS_FATAL("No suitable TESFEMReactionAdaptor found. Aborting.");
    return nullptr;
}

TESFEMReactionAdaptorAdsorption::TESFEMReactionAdaptorAdsorption(
    TESLocalAssemblerData const& data)
    // caution fragile: this relies in this constructor b eing called __after__
    // data.solid_density has been properly set up!
    : bounds_violation_(data.solid_density.size(), false),
      d_(data)
{
    assert(dynamic_cast<Adsorption::AdsorptionReaction const*>(
               data.ap.react_sys.get()) != nullptr &&
           "Reactive system has wrong type.");
    assert(!bounds_violation_.empty());
}

ReactionRate
TESFEMReactionAdaptorAdsorption::initReaction_slowDownUndershootStrategy(
    const unsigned int_pt)
{
    assert(d_.ap.number_of_try_of_iteration <= 20);

    const double loading = Adsorption::AdsorptionReaction::getLoading(
        d_.solid_density_prev_ts[int_pt], d_.ap.rho_SR_dry);

    double react_rate_R =
        d_.ap.react_sys->getReactionRate(d_.p_V, d_.T, d_.ap.M_react, loading) *
        d_.ap.rho_SR_dry;

    // set reaction rate based on current damping factor
    react_rate_R = (reaction_damping_factor_ > 1e-3)
                       ? reaction_damping_factor_ * react_rate_R
                       : 0.0;

    if (d_.p_V <
            0.01 * Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(
                       d_.T) &&
        react_rate_R > 0.0)
    {
        react_rate_R = 0.0;
    }
    else if (d_.p_V < 100.0 ||
             d_.p_V < 0.05 * Adsorption::AdsorptionReaction::
                                 getEquilibriumVapourPressure(d_.T))
    {
        // use equilibrium reaction for dry regime

        // in the case of zeroth try in zeroth iteration: p_V_ and loading are
        // the values
        // at the end of the previous timestep

        const double pV_eq = estimateAdsorptionEquilibrium(d_.p_V, loading);
        // TODO [CL]: it would be more correct to subtract pV from the previous
        // timestep here
        const double delta_pV = pV_eq - d_.p_V;
        const double delta_rhoV =
            delta_pV * d_.ap.M_react /
            MaterialLib::PhysicalConstant::IdealGasConstant / d_.T *
            d_.ap.poro;
        const double delta_rhoSR = delta_rhoV / (d_.ap.poro - 1.0);
        double react_rate_R2 = delta_rhoSR / d_.ap.delta_t;

        if (bounds_violation_[int_pt])
        {
            react_rate_R2 *= 0.5;
        }

        // 0th try: make sure reaction is not slower than allowed by local
        // estimation
        // nth try: make sure reaction is not made faster by local estimation
        if ((d_.ap.number_of_try_of_iteration == 1 &&
             std::abs(react_rate_R2) > std::abs(react_rate_R)) ||
            (d_.ap.number_of_try_of_iteration > 1 &&
             std::abs(react_rate_R2) < std::abs(react_rate_R)))
        {
            react_rate_R = react_rate_R2;
        }
    }

    // smooth out readjustment of reaction rate
    if (d_.ap.iteration_in_current_timestep > 4)
    {
        if (d_.ap.iteration_in_current_timestep <= 9)
        {
            // update reaction rate for for five iterations
            const auto N = d_.ap.iteration_in_current_timestep - 4;

            // take average s.t. does not oscillate so much
            react_rate_R = 1.0 / (1.0 + N) *
                           (N * d_.reaction_rate[int_pt] + 1.0 * react_rate_R);
        }
        else
        {
            // afterwards no update anymore
            react_rate_R = d_.reaction_rate[int_pt];
        }
    }

    if (d_.ap.number_of_try_of_iteration > 1)
    {
        // assert that within tries reaction does not get faster
        // (e.g. due to switch equilibrium reaction <--> kinetic reaction)

        // factor of 0.9*N: in fact, even slow down reaction over tries
        const double r = std::pow(0.9, d_.ap.number_of_try_of_iteration - 1) *
                         d_.reaction_rate[int_pt];
        if (std::abs(react_rate_R) > std::abs(r))
        {
            react_rate_R = r;
        }
    }

    return {react_rate_R,
            d_.solid_density_prev_ts[int_pt] + react_rate_R * d_.ap.delta_t};
}

double TESFEMReactionAdaptorAdsorption::estimateAdsorptionEquilibrium(
    const double p_V0, const double C0) const
{
    auto f = [this, p_V0, C0](double pV) -> double {
        // pV0 := p_V_
        const double C_eq =
            d_.ap.react_sys->getEquilibriumLoading(pV, d_.T, d_.ap.M_react);
        return (pV - p_V0) * d_.ap.M_react /
                   MaterialLib::PhysicalConstant::IdealGasConstant / d_.T *
                   d_.ap.poro +
               (1.0 - d_.ap.poro) * (C_eq - C0) * d_.ap.rho_SR_dry;
    };

    // range where to search for roots of f
    const double C_eq0 =
        d_.ap.react_sys->getEquilibriumLoading(p_V0, d_.T, d_.ap.M_react);
    const double limit =
        (C_eq0 > C0)
            ? 1e-8
            : Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(
                  d_.T);

    // search for roots
    auto rf = MathLib::Nonlinear::makeRegulaFalsi<MathLib::Nonlinear::Pegasus>(
        f, p_V0, limit);
    rf.step(3);

    // set vapour pressure
    return rf.getResult();
}

bool TESFEMReactionAdaptorAdsorption::checkBounds(
    std::vector<double> const& local_x,
    std::vector<double> const& local_x_prev_ts)
{
    double alpha = 1.0;

    const double min_xmV = 1e-6;
    const std::size_t nnodes = local_x.size() / NODAL_DOF;
    const std::size_t xmV_offset = COMPONENT_ID_MASS_FRACTION * nnodes;

    for (std::size_t i = 0; i < nnodes; ++i)
    {
        auto const xnew = local_x[xmV_offset + i];
        auto const xold = local_x_prev_ts[xmV_offset + i];

        if (xnew < min_xmV)
        {
            const auto a = xold / (xold - xnew);
            alpha = std::min(alpha, a);
            bounds_violation_[i] = true;
        }
        else if (xnew > 1.0)
        {
            const auto a = xold / (xnew - xold);
            alpha = std::min(alpha, a);
            bounds_violation_[i] = true;
        }
        else
        {
            bounds_violation_[i] = false;
        }
    }

    assert(alpha > 0.0);

    if (alpha != 1.0)
    {
        if (alpha > 0.5)
        {
            alpha = 0.5;
        }
        if (alpha < 0.05)
        {
            alpha = 0.05;
        }

        if (d_.ap.number_of_try_of_iteration <= 3)
        {
            reaction_damping_factor_ *= sqrt(alpha);
        }
        else
        {
            reaction_damping_factor_ *= alpha;
        }
    }

    return alpha == 1.0;
}

void TESFEMReactionAdaptorAdsorption::preZerothTryAssemble()
{
    if (reaction_damping_factor_ < 1e-3)
    {
        reaction_damping_factor_ = 1e-3;
    }

    reaction_damping_factor_ = std::min(std::sqrt(reaction_damping_factor_),
                                        10.0 * reaction_damping_factor_);
}

TESFEMReactionAdaptorInert::TESFEMReactionAdaptorInert(
    TESLocalAssemblerData const& data)
    : d_(data)
{
}

ReactionRate TESFEMReactionAdaptorInert::initReaction(const unsigned int_pt)
{
    return {0.0, d_.solid_density_prev_ts[int_pt]};
    // d_.qR_ = 0.0;
    // d_.reaction_rate_[int_pt] = 0.0;
}

TESFEMReactionAdaptorSinusoidal::TESFEMReactionAdaptorSinusoidal(
    TESLocalAssemblerData const& data)
    : d_(data)
{
    assert(dynamic_cast<Adsorption::ReactionSinusoidal const*>(
               data.ap.react_sys.get()) != nullptr &&
           "Reactive system has wrong type.");
}

ReactionRate TESFEMReactionAdaptorSinusoidal::initReaction(
    const unsigned /*int_pt*/)
{
    const double t = d_.ap.current_time;

    // Cf. OGS5
    const double rhoSR0 = 1.0;
    const double rhoTil = 0.1;
    const double omega = 2.0 * 3.1416;
    const double poro = d_.ap.poro;

    return {rhoTil * omega * cos(omega * t) / (1.0 - poro),
            rhoSR0 + rhoTil * std::sin(omega * t) / (1.0 - poro)};
}

TESFEMReactionAdaptorCaOH2::TESFEMReactionAdaptorCaOH2(
    TESLocalAssemblerData const& data)
    : d_(data),
      react_(dynamic_cast<Adsorption::ReactionCaOH2&>(*data.ap.react_sys))
{
    ode_solver_ = MathLib::ODE::createODESolver<1>(react_.getOdeSolverConfig());
    // TODO invalidate config

    ode_solver_->setTolerance(1e-10, 1e-10);

    auto f = [this](const double /*t*/,
                    MathLib::ODE::MappedConstVector<1> const y,
                    MathLib::ODE::MappedVector<1>
                        ydot) -> bool {
        ydot[0] = react_.getReactionRate(y[0]);
        return true;
    };

    ode_solver_->setFunction(f, nullptr);
}

ReactionRate TESFEMReactionAdaptorCaOH2::initReaction(const unsigned int int_pt)
{
    // TODO if the first holds, the second also has to hold
    if (d_.ap.iteration_in_current_timestep > 1 ||
        d_.ap.number_of_try_of_iteration > 1)
    {
        return {d_.reaction_rate[int_pt], d_.solid_density[int_pt]};
    }

    // TODO: double check!
    // const double xv_NR  = SolidProp->non_reactive_solid_volume_fraction;
    // const double rho_NR = SolidProp->non_reactive_solid_density;
    const double xv_NR = 0.0;
    const double rho_NR = 0.0;

    const double t0 = 0.0;
    const double y0 =
        (d_.solid_density_prev_ts[int_pt] - xv_NR * rho_NR) / (1.0 - xv_NR);

    const double t_end = d_.ap.delta_t;

    react_.updateParam(d_.T, d_.p, d_.vapour_mass_fraction,
                       d_.solid_density_prev_ts[int_pt]);

    ode_solver_->setIC(t0, {y0});
    ode_solver_->preSolve();
    ode_solver_->solve(t_end);

    const double time_reached = ode_solver_->getTime();
    (void)time_reached;
    assert(std::abs(t_end - time_reached) <
           std::numeric_limits<double>::epsilon());

    auto const& y_new = ode_solver_->getSolution();
    auto const& y_dot_new = ode_solver_->getYDot(t_end, y_new);

    double rho_react;

    // cut off when limits are reached
    if (y_new[0] < ProcessLib::TES::TESFEMReactionAdaptorCaOH2::React::rho_low)
    {
        rho_react = ProcessLib::TES::TESFEMReactionAdaptorCaOH2::React::rho_low;
    }
    else if (y_new[0] >
             ProcessLib::TES::TESFEMReactionAdaptorCaOH2::React::rho_up)
    {
        rho_react = ProcessLib::TES::TESFEMReactionAdaptorCaOH2::React::rho_up;
    }
    else
    {
        rho_react = y_new[0];
    }

    return {y_dot_new[0] * (1.0 - xv_NR),
            (1.0 - xv_NR) * rho_react + xv_NR * rho_NR};
}

}  // namespace TES
}  // namespace ProcessLib
