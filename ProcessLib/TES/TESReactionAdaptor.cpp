/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cassert>

#include <logog/include/logog.hpp>

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
    : _bounds_violation(data.solid_density.size(), false),
      _d(data)
{
    assert(dynamic_cast<Adsorption::AdsorptionReaction const*>(
               data.ap.react_sys.get()) != nullptr &&
           "Reactive system has wrong type.");
    assert(!_bounds_violation.empty());
}

ReactionRate
TESFEMReactionAdaptorAdsorption::initReaction_slowDownUndershootStrategy(
    const unsigned int_pt)
{
    assert(_d.ap.number_of_try_of_iteration <= 20);

    const double loading = Adsorption::AdsorptionReaction::getLoading(
        _d.solid_density_prev_ts[int_pt], _d.ap.rho_SR_dry);

    double react_rate_R =
        _d.ap.react_sys->getReactionRate(_d.p_V, _d.T, _d.ap.M_react, loading) *
        _d.ap.rho_SR_dry;

    // set reaction rate based on current damping factor
    react_rate_R = (_reaction_damping_factor > 1e-3)
                       ? _reaction_damping_factor * react_rate_R
                       : 0.0;

    if (_d.p_V <
            0.01 * Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(
                       _d.T) &&
        react_rate_R > 0.0)
    {
        react_rate_R = 0.0;
    }
    else if (_d.p_V < 100.0 ||
             _d.p_V < 0.05 * Adsorption::AdsorptionReaction::
                                 getEquilibriumVapourPressure(_d.T))
    {
        // use equilibrium reaction for dry regime

        // in the case of zeroth try in zeroth iteration: _p_V and loading are
        // the values
        // at the end of the previous timestep

        const double pV_eq = estimateAdsorptionEquilibrium(_d.p_V, loading);
        // TODO [CL]: it would be more correct to subtract pV from the previous
        // timestep here
        const double delta_pV = pV_eq - _d.p_V;
        const double delta_rhoV =
            delta_pV * _d.ap.M_react /
            MaterialLib::PhysicalConstant::IdealGasConstant / _d.T *
            _d.ap.poro;
        const double delta_rhoSR = delta_rhoV / (_d.ap.poro - 1.0);
        double react_rate_R2 = delta_rhoSR / _d.ap.delta_t;

        if (_bounds_violation[int_pt])
        {
            react_rate_R2 *= 0.5;
        }

        // 0th try: make sure reaction is not slower than allowed by local
        // estimation
        // nth try: make sure reaction is not made faster by local estimation
        if ((_d.ap.number_of_try_of_iteration == 1 &&
             std::abs(react_rate_R2) > std::abs(react_rate_R)) ||
            (_d.ap.number_of_try_of_iteration > 1 &&
             std::abs(react_rate_R2) < std::abs(react_rate_R)))
        {
            react_rate_R = react_rate_R2;
        }
    }

    // smooth out readjustment of reaction rate
    if (_d.ap.iteration_in_current_timestep > 4)
    {
        if (_d.ap.iteration_in_current_timestep <= 9)
        {
            // update reaction rate for for five iterations
            const auto N = _d.ap.iteration_in_current_timestep - 4;

            // take average s.t. does not oscillate so much
            react_rate_R = 1.0 / (1.0 + N) *
                           (N * _d.reaction_rate[int_pt] + 1.0 * react_rate_R);
        }
        else
        {
            // afterwards no update anymore
            react_rate_R = _d.reaction_rate[int_pt];
        }
    }

    if (_d.ap.number_of_try_of_iteration > 1)
    {
        // assert that within tries reaction does not get faster
        // (e.g. due to switch equilibrium reaction <--> kinetic reaction)

        // factor of 0.9*N: in fact, even slow down reaction over tries
        const double r = std::pow(0.9, _d.ap.number_of_try_of_iteration - 1) *
                         _d.reaction_rate[int_pt];
        if (std::abs(react_rate_R) > std::abs(r))
        {
            react_rate_R = r;
        }
    }

    return {react_rate_R,
            _d.solid_density_prev_ts[int_pt] + react_rate_R * _d.ap.delta_t};
}

double TESFEMReactionAdaptorAdsorption::estimateAdsorptionEquilibrium(
    const double p_V0, const double C0) const
{
    auto f = [this, p_V0, C0](double pV) -> double {
        // pV0 := _p_V
        const double C_eq =
            _d.ap.react_sys->getEquilibriumLoading(pV, _d.T, _d.ap.M_react);
        return (pV - p_V0) * _d.ap.M_react /
                   MaterialLib::PhysicalConstant::IdealGasConstant / _d.T *
                   _d.ap.poro +
               (1.0 - _d.ap.poro) * (C_eq - C0) * _d.ap.rho_SR_dry;
    };

    // range where to search for roots of f
    const double C_eq0 =
        _d.ap.react_sys->getEquilibriumLoading(p_V0, _d.T, _d.ap.M_react);
    const double limit =
        (C_eq0 > C0)
            ? 1e-8
            : Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(
                  _d.T);

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
            _bounds_violation[i] = true;
        }
        else if (xnew > 1.0)
        {
            const auto a = xold / (xnew - xold);
            alpha = std::min(alpha, a);
            _bounds_violation[i] = true;
        }
        else
        {
            _bounds_violation[i] = false;
        }
    }

    assert(alpha > 0.0);

    if (alpha != 1.0)
    {
        if (alpha > 0.5)
            alpha = 0.5;
        if (alpha < 0.05)
            alpha = 0.05;

        if (_d.ap.number_of_try_of_iteration <= 3)
        {
            _reaction_damping_factor *= sqrt(alpha);
        }
        else
        {
            _reaction_damping_factor *= alpha;
        }
    }

    return alpha == 1.0;
}

void TESFEMReactionAdaptorAdsorption::preZerothTryAssemble()
{
    if (_reaction_damping_factor < 1e-3)
        _reaction_damping_factor = 1e-3;

    _reaction_damping_factor = std::min(std::sqrt(_reaction_damping_factor),
                                        10.0 * _reaction_damping_factor);
}

TESFEMReactionAdaptorInert::TESFEMReactionAdaptorInert(
    TESLocalAssemblerData const& data)
    : _d(data)
{
}

ReactionRate TESFEMReactionAdaptorInert::initReaction(const unsigned int_pt)
{
    return {0.0, _d.solid_density_prev_ts[int_pt]};
    // _d._qR = 0.0;
    // _d._reaction_rate[int_pt] = 0.0;
}

TESFEMReactionAdaptorSinusoidal::TESFEMReactionAdaptorSinusoidal(
    TESLocalAssemblerData const& data)
    : _d(data)
{
    assert(dynamic_cast<Adsorption::ReactionSinusoidal const*>(
               data.ap.react_sys.get()) != nullptr &&
           "Reactive system has wrong type.");
}

ReactionRate TESFEMReactionAdaptorSinusoidal::initReaction(
    const unsigned /*int_pt*/)
{
    const double t = _d.ap.current_time;

    // Cf. OGS5
    const double rhoSR0 = 1.0;
    const double rhoTil = 0.1;
    const double omega = 2.0 * 3.1416;
    const double poro = _d.ap.poro;

    return {rhoTil * omega * cos(omega * t) / (1.0 - poro),
            rhoSR0 + rhoTil * std::sin(omega * t) / (1.0 - poro)};
}

TESFEMReactionAdaptorCaOH2::TESFEMReactionAdaptorCaOH2(
    TESLocalAssemblerData const& data)
    : _d(data),
      _react(dynamic_cast<Adsorption::ReactionCaOH2&>(*data.ap.react_sys))
{
    _ode_solver = MathLib::ODE::createODESolver<1>(_react.getOdeSolverConfig());
    // TODO invalidate config

    _ode_solver->setTolerance(1e-10, 1e-10);

    auto f = [this](const double /*t*/,
                    MathLib::ODE::MappedConstVector<1> const y,
                    MathLib::ODE::MappedVector<1>
                        ydot) -> bool {
        ydot[0] = _react.getReactionRate(y[0]);
        return true;
    };

    _ode_solver->setFunction(f, nullptr);
}

ReactionRate TESFEMReactionAdaptorCaOH2::initReaction(const unsigned int int_pt)
{
    // TODO if the first holds, the second also has to hold
    if (_d.ap.iteration_in_current_timestep > 1 ||
        _d.ap.number_of_try_of_iteration > 1)
    {
        return {_d.reaction_rate[int_pt], _d.solid_density[int_pt]};
    }

    // TODO: double check!
    // const double xv_NR  = SolidProp->non_reactive_solid_volume_fraction;
    // const double rho_NR = SolidProp->non_reactive_solid_density;
    const double xv_NR = 0.0;
    const double rho_NR = 0.0;

    const double t0 = 0.0;
    const double y0 =
        (_d.solid_density_prev_ts[int_pt] - xv_NR * rho_NR) / (1.0 - xv_NR);

    const double t_end = _d.ap.delta_t;

    _react.updateParam(_d.T, _d.p, _d.vapour_mass_fraction,
                       _d.solid_density_prev_ts[int_pt]);

    _ode_solver->setIC(t0, {y0});
    _ode_solver->preSolve();
    _ode_solver->solve(t_end);

    const double time_reached = _ode_solver->getTime();
    (void)time_reached;
    assert(std::abs(t_end - time_reached) <
           std::numeric_limits<double>::epsilon());

    auto const& y_new = _ode_solver->getSolution();
    auto const& y_dot_new = _ode_solver->getYDot(t_end, y_new);

    double rho_react;

    // cut off when limits are reached
    if (y_new[0] < _react.rho_low)
        rho_react = _react.rho_low;
    else if (y_new[0] > _react.rho_up)
        rho_react = _react.rho_up;
    else
        rho_react = y_new[0];

    return {y_dot_new[0] * (1.0 - xv_NR),
            (1.0 - xv_NR) * rho_react + xv_NR * rho_NR};
}

}  // namespace TES
}  // namespace ProcessLib
