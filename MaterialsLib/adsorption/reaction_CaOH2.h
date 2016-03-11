/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ArrayRef.h"

#include "reaction.h"
#include "adsorption.h"


namespace ProcessLib
{
template<typename>
class TESFEMReactionAdaptorCaOH2;
}

namespace Ads
{

class ReactionCaOH2 final : public Reaction
{
public:
    explicit ReactionCaOH2(BaseLib::ConfigTreeNew const& conf)
        : ode_solver_config{conf.getConfSubtree("ode_solver_config")}
    {
        /*auto const param = conf.get_optional<double>("reaction_enthalpy");
        if (param) {
            _enthalpy = *param;
        } else {
            ERR("<reaction_enthalpy> not specified.");
            std::abort();
        }*/
    }

    double get_enthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
                        const double /*M_Ads*/) const override;

    double get_reaction_rate(const double /*p_Ads*/, const double /*T_Ads*/, const double /*M_Ads*/,
                             const double /*loading*/) const override;

    const BaseLib::ConfigTreeNew& getOdeSolverConfig() const { return ode_solver_config; }


    void eval(double /*t*/,
              BaseLib::ArrayRef<const double, 1> const& y,
              BaseLib::ArrayRef<double, 1>& dydx);
    void update_param(double T_solid,
                      double p_gas,
                      double x_react,
                      double rho_s_initial);

private:
    void calculate_qR();
    void set_chemical_equilibrium();
    double Ca_hydration();

    static constexpr double R = Ads::GAS_CONST;  // [J/mol/K]
    double rho_s;                 // solid phase density
    double p_gas;                 // gas phase pressure in unit bar
    double p_r_g;                 // pressure of H2O on gas phase;
    double p_eq       = 1.0;      // equilibrium pressure; // unit in bar
    double T_eq;                  // equilibrium temperature;
    double T_s;                   // solid phase temperature;
    double qR;                    // rate of solid density change;
    double x_react;               // mass fraction of water in gas phase;
    double X_D;                   // mass fraction of dehydration (CaO) in the solid phase;
    double X_H;                   // mass fraction of hydration in the solid phase;
    static constexpr double reaction_enthalpy = -1.12e+05; // in J/mol; negative for exothermic composition reaction
    static constexpr double reaction_entropy  = -143.5;    // in J/mol/K
    static constexpr double M_carrier = Ads::M_N2;         // inert component molar mass
    static constexpr double M_react   = Ads::M_H2O;        // reactive component molar mass

    static constexpr double tol_l   = 1e-4;
    static constexpr double tol_u   = 1.0 - 1e-4;
    static constexpr double tol_rho = 0.1;

    const BaseLib::ConfigTreeNew ode_solver_config;

    template<typename>
    friend class ProcessLib::TESFEMReactionAdaptorCaOH2;

public:
    static constexpr double rho_low = 1656.0;              // lower density limit
    static constexpr double rho_up = 2200.0;               // upper density limit
};

}
