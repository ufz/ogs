/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <memory>
#include "TESFEM-data.h"

#include "MaterialsLib/adsorption/reaction_CaOH2.h"
#include "MathLib/ODE/OdeSolverFactory.h"



namespace ProcessLib
{
namespace TES
{

template<typename Traits>
class LADataNoTpl;

struct ReactionRate
{
    const double reaction_rate;
    const double solid_density;
};

template<typename Traits>
class TESFEMReactionAdaptor
{
public:
    virtual bool
    checkBounds(std::vector<double> const& localX,
                std::vector<double> const& localX_pts) = 0;

    virtual ReactionRate
    initReaction(const unsigned int_pt) = 0;

    virtual void
    preZerothTryAssemble() = 0;

    // TODO: remove
    virtual double getReactionDampingFactor() const {
        return -1.0;
    }

    virtual ~TESFEMReactionAdaptor() = default;

    static std::unique_ptr<TESFEMReactionAdaptor<Traits> >
    newInstance(LADataNoTpl<Traits>& data);
};


template<typename Traits>
class TESFEMReactionAdaptorAdsorption final : public TESFEMReactionAdaptor<Traits>
{
public:
    explicit TESFEMReactionAdaptorAdsorption(LADataNoTpl<Traits>& data);

    bool checkBounds(const std::vector<double> &localX, const std::vector<double> &localX_pts)
    override;

    ReactionRate initReaction(const unsigned int_pt) override {
        return initReaction_slowDownUndershootStrategy(int_pt);
    }

    void preZerothTryAssemble() override;

    // TODO: get rid of
    double getReactionDampingFactor() const override {
        return _reaction_damping_factor;
    }

private:
    ReactionRate initReaction_slowDownUndershootStrategy(const unsigned int_pt);

    /// returns estimated equilibrium vapour pressure
    /// based on a local (i.e. no diffusion/advection) balance
    /// of adsorbate loading and vapour partial pressure
    double estimateAdsorptionEquilibrium(const double p_V0, const double C0) const;

    double _reaction_damping_factor = 1.0;
    std::vector<bool> _bounds_violation;

    LADataNoTpl<Traits>& _data;
};


template<typename Traits>
class TESFEMReactionAdaptorInert final : public TESFEMReactionAdaptor<Traits>
{
public:
    explicit TESFEMReactionAdaptorInert(LADataNoTpl<Traits>&);

    bool checkBounds(const std::vector<double> &, const std::vector<double> &)
    override {
        return true;
    }

    ReactionRate initReaction(const unsigned int_pt) override
    {
        return { 0.0, _data._solid_density_prev_ts[int_pt] };
        // _data._qR = 0.0;
        // _data._reaction_rate[int_pt] = 0.0;
    }

    void preZerothTryAssemble() override
    {}

private:
    LADataNoTpl<Traits>& _data;
};


template<typename Traits>
class TESFEMReactionAdaptorSinusoidal final : public TESFEMReactionAdaptor<Traits>
{
public:
    explicit TESFEMReactionAdaptorSinusoidal(LADataNoTpl<Traits>& data);

    bool checkBounds(const std::vector<double> &, const std::vector<double> &)
    override {
        return true;
    }

    ReactionRate initReaction(const unsigned) override;

    void preZerothTryAssemble() override
    {}

private:
    LADataNoTpl<Traits>& _data;
};


template<typename Traits>
class TESFEMReactionAdaptorCaOH2 final : public TESFEMReactionAdaptor<Traits>
{
public:
    explicit TESFEMReactionAdaptorCaOH2(LADataNoTpl<Traits>& data);

    bool checkBounds(const std::vector<double> &, const std::vector<double> &)
    override {
        return true;
    }

    ReactionRate initReaction(const unsigned) override;

    void preZerothTryAssemble() override
    {}

private:
    using Data = LADataNoTpl<Traits>;
    using React = Ads::ReactionCaOH2;
    Data& _data;
    React& _react;

    std::unique_ptr<MathLib::OdeSolver<1, React> > _ode_solver;

    static bool odeRhs(const double /*t*/,
                       BaseLib::ArrayRef<const double, 1> y,
                       BaseLib::ArrayRef<double, 1> ydot,
                       React& reaction);
};

}
}

#include "TESFEMReactionAdaptor-impl.h"
