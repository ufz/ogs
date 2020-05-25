/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "materiallib_export.h"

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
        ode_solver_config_{conf.getConfigSubtree("ode_solver_config")}
    {}

    double getEnthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
                        const double /*M_Ads*/) const override;

    double getReactionRate(const double /*p_Ads*/, const double /*T_Ads*/, const double /*M_Ads*/,
                             const double /*loading*/) const override;

    const BaseLib::ConfigTree& getOdeSolverConfig() const { return ode_solver_config_; }

    // TODO merge with getReactionRate() above
    double getReactionRate(double const solid_density);

    void updateParam(double T_solid,
                      double p_gas_,
                      double x_react_,
                      double rho_s_initial);

private:
    void calculateQR();
    void setChemicalEquilibrium();
    double CaHydration();

    double rho_s_;           //!< solid phase density
    double p_gas_;           //!< gas phase pressure in unit bar
    double p_r_g_;           //!< pressure of H2O on gas phase
    double p_eq_ = 1.0;      //!< equilibrium pressure in bar
    double T_eq_;            //!< equilibrium temperature
    double T_s_;             //!< solid phase temperature
    double qR_;              //!< rate of solid density change
    double x_react_;         //!< mass fraction of water in gas phase
    double X_D_;             //!< mass fraction of dehydration (CaO) in the solid phase
    double X_H_;             //!< mass fraction of hydration in the solid phase

    //! reaction enthalpy in J/mol; negative for exothermic composition reaction
    static const double reaction_enthalpy_;
    static const double reaction_entropy_; //!< reaction entropy in J/mol/K
    static const double M_carrier_;        //!< inert component molar mass
    static const double M_react_;          //!< reactive component molar mass

    static const double tol_l_;
    static const double tol_u_;
    static const double tol_rho_;

    const BaseLib::ConfigTree ode_solver_config_;

    template<typename>
    friend class ProcessLib::TESFEMReactionAdaptorCaOH2;

public:
    static MATERIALLIB_EXPORT const double rho_low; //! lower density limit
    static MATERIALLIB_EXPORT const double rho_up;  //! upper density limit
};

}  // namespace Adsorption
