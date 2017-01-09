/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_TES_TESREACTIONADAPTOR_H
#define PROCESSLIB_TES_TESREACTIONADAPTOR_H

#include <memory>
#include <vector>

#include "MaterialLib/Adsorption/ReactionCaOH2.h"
#include "MathLib/ODE/ODESolver.h"

namespace ProcessLib
{
namespace TES
{
struct TESLocalAssemblerData;

struct ReactionRate
{
    const double reaction_rate;
    const double solid_density;
};

class TESFEMReactionAdaptor
{
public:
    virtual bool checkBounds(std::vector<double> const& /*local_x*/,
                             std::vector<double> const& /*local_x_prev_ts*/)
    {
        return true;  // by default accept everything
    }

    virtual ReactionRate initReaction(const unsigned int_pt) = 0;

    virtual void preZerothTryAssemble() {}
    // TODO: remove
    virtual double getReactionDampingFactor() const { return -1.0; }
    virtual ~TESFEMReactionAdaptor() = default;

    static std::unique_ptr<TESFEMReactionAdaptor> newInstance(
        TESLocalAssemblerData const& data);
};

class TESFEMReactionAdaptorAdsorption final : public TESFEMReactionAdaptor
{
public:
    explicit TESFEMReactionAdaptorAdsorption(TESLocalAssemblerData const& data);

    bool checkBounds(std::vector<double> const& local_x,
                     std::vector<double> const& local_x_prev_ts) override;

    ReactionRate initReaction(const unsigned int_pt) override
    {
        return initReaction_slowDownUndershootStrategy(int_pt);
    }

    void preZerothTryAssemble() override;

    // TODO: get rid of
    double getReactionDampingFactor() const override
    {
        return _reaction_damping_factor;
    }

private:
    ReactionRate initReaction_slowDownUndershootStrategy(const unsigned int_pt);

    /// returns estimated equilibrium vapour pressure
    /// based on a local (i.e. no diffusion/advection) balance
    /// of adsorbate loading and vapour partial pressure
    double estimateAdsorptionEquilibrium(const double p_V0,
                                         const double C0) const;

    double _reaction_damping_factor = 1.0;
    std::vector<bool> _bounds_violation;

    TESLocalAssemblerData const& _d;
};

class TESFEMReactionAdaptorInert final : public TESFEMReactionAdaptor
{
public:
    explicit TESFEMReactionAdaptorInert(TESLocalAssemblerData const&);

    ReactionRate initReaction(const unsigned int_pt) override;

private:
    TESLocalAssemblerData const& _d;
};

class TESFEMReactionAdaptorSinusoidal final : public TESFEMReactionAdaptor
{
public:
    explicit TESFEMReactionAdaptorSinusoidal(TESLocalAssemblerData const& data);

    ReactionRate initReaction(const unsigned) override;

private:
    TESLocalAssemblerData const& _d;
};

class TESFEMReactionAdaptorCaOH2 final : public TESFEMReactionAdaptor
{
public:
    explicit TESFEMReactionAdaptorCaOH2(TESLocalAssemblerData const& data);

    ReactionRate initReaction(const unsigned) override;

private:
    using Data = TESLocalAssemblerData;
    using React = Adsorption::ReactionCaOH2;
    Data const& _d;
    React& _react;

    std::unique_ptr<MathLib::ODE::ODESolver<1>> _ode_solver;

    static bool odeRhs(const double /*t*/,
                       MathLib::ODE::MappedConstVector<1> const y,
                       MathLib::ODE::MappedVector<1>
                           ydot,
                       React& reaction);
};

}  // namespace TES
}  // namespace ProcessLib
#endif  // PROCESSLIB_TES_TESREACTIONADAPTOR_H
