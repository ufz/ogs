/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "Reaction.h"
#include "Adsorption.h"

namespace ProcessLib
{
template<typename>
class TESFEMReactionAdaptorCaOH2;
}

namespace Adsorption
{

class ReactionCaOH2 final : public Reaction
{
public:
    explicit ReactionCaOH2(BaseLib::ConfigTree const& conf)
        : ode_solver_config{conf.getConfSubtree("ode_solver_config")}
    {}

    double getEnthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
                        const double /*M_Ads*/) const override;

    double getReactionRate(const double /*p_Ads*/, const double /*T_Ads*/, const double /*M_Ads*/,
                             const double /*loading*/) const override;

    const BaseLib::ConfigTree& getOdeSolverConfig() const { return ode_solver_config; }

    // TODO merge with getReactionRate() above
    double getReactionRate(double const solid_density);

    void updateParam(double T_solid,
                      double p_gas,
                      double x_react,
                      double rho_s_initial);

private:
    void calculateQR();
    void setChemicalEquilibrium();
    double CaHydration();

    static const double R;  // [J/mol/K]
    double rho_s;           // solid phase density
    double p_gas;           // gas phase pressure in unit bar
    double p_r_g;           // pressure of H2O on gas phase;
    double p_eq = 1.0;      // equilibrium pressure; // unit in bar
    double T_eq;            // equilibrium temperature;
    double T_s;             // solid phase temperature;
    double qR;              // rate of solid density change;
    double x_react;         // mass fraction of water in gas phase;
    double X_D;             // mass fraction of dehydration (CaO) in the solid phase;
    double X_H;             // mass fraction of hydration in the solid phase;
    static const double reaction_enthalpy; // in J/mol; negative for exothermic composition reaction
    static const double reaction_entropy; // in J/mol/K
    static const double M_carrier;        // inert component molar mass
    static const double M_react;          // reactive component molar mass

    static const double tol_l;
    static const double tol_u;
    static const double tol_rho;

    const BaseLib::ConfigTree ode_solver_config;

    template<typename>
    friend class ProcessLib::TESFEMReactionAdaptorCaOH2;

public:
    static const double rho_low; // lower density limit
    static const double rho_up;  // upper density limit
};

}
