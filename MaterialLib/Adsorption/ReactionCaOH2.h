/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    explicit ReactionCaOH2(BaseLib::ConfigTree const& conf) :
        //! \ogs_file_param{material__adsorption__reaction__CaOH2__ode_solver_config}
        _ode_solver_config{conf.getConfigSubtree("ode_solver_config")}
    {}

    double getEnthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
                        const double /*M_Ads*/) const override;

    double getReactionRate(const double /*p_Ads*/, const double /*T_Ads*/, const double /*M_Ads*/,
                             const double /*loading*/) const override;

    const BaseLib::ConfigTree& getOdeSolverConfig() const { return _ode_solver_config; }

    // TODO merge with getReactionRate() above
    double getReactionRate(double const solid_density);

    void updateParam(double T_solid,
                      double _p_gas,
                      double _x_react,
                      double rho_s_initial);

private:
    void calculateQR();
    void setChemicalEquilibrium();
    double CaHydration();

    double _rho_s;           //!< solid phase density
    double _p_gas;           //!< gas phase pressure in unit bar
    double _p_r_g;           //!< pressure of H2O on gas phase
    double _p_eq = 1.0;      //!< equilibrium pressure in bar
    double _T_eq;            //!< equilibrium temperature
    double _T_s;             //!< solid phase temperature
    double _qR;              //!< rate of solid density change
    double _x_react;         //!< mass fraction of water in gas phase
    double _X_D;             //!< mass fraction of dehydration (CaO) in the solid phase
    double _X_H;             //!< mass fraction of hydration in the solid phase

    //! reaction enthalpy in J/mol; negative for exothermic composition reaction
    static const double _reaction_enthalpy;
    static const double _reaction_entropy; //!< reaction entropy in J/mol/K
    static const double _M_carrier;        //!< inert component molar mass
    static const double _M_react;          //!< reactive component molar mass

    static const double _tol_l;
    static const double _tol_u;
    static const double _tol_rho;

    const BaseLib::ConfigTree _ode_solver_config;

    template<typename>
    friend class ProcessLib::TESFEMReactionAdaptorCaOH2;

public:
    static const double rho_low; //! lower density limit
    static const double rho_up;  //! upper density limit
};

}
