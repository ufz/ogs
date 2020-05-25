/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ReactionCaOH2.h"
#include <cassert>
#include "Adsorption.h"
#include "BaseLib/Error.h"
#include "MaterialLib/PhysicalConstant.h"

namespace Adsorption
{

const double ReactionCaOH2::reaction_enthalpy_ = -1.12e+05;
const double ReactionCaOH2::reaction_entropy_  = -143.5;
const double ReactionCaOH2::M_carrier_ =
    MaterialLib::PhysicalConstant::MolarMass::N2;
const double ReactionCaOH2::M_react_ =
    MaterialLib::PhysicalConstant::MolarMass::Water;

const double ReactionCaOH2::tol_l_   = 1e-4;
const double ReactionCaOH2::tol_u_   = 1.0 - 1e-4;
const double ReactionCaOH2::tol_rho_ = 0.1;

const double ReactionCaOH2::rho_low = 1656.0;
const double ReactionCaOH2::rho_up = 2200.0;

double ReactionCaOH2::getEnthalpy(const double /*p_Ads*/,
                                  const double /*T_Ads*/,
                                  const double /*M_Ads*/) const
{
    return - reaction_enthalpy_/M_react_;
}

double ReactionCaOH2::getReactionRate(const double /*p_Ads*/,
                                      const double /*T_Ads*/,
                                      const double /*M_Ads*/,
                                      const double /*loading*/) const
{
    OGS_FATAL("get_reaction_rate do not call directly");
}


double ReactionCaOH2::getReactionRate(double const solid_density)
{
    rho_s_ = solid_density;
    calculateQR();
    return qR_;
}

void ReactionCaOH2::updateParam(
        double T_solid,
        double p_gas,
        double x_react,
        double rho_s_initial)
{
    T_s_     = T_solid;
    p_gas_ = p_gas / 1e5; // convert Pa to bar
    x_react_ = x_react;
    rho_s_   = rho_s_initial;
}

void ReactionCaOH2::calculateQR()
{
    // Convert mass fraction into mole fraction
    const double mol_frac_react = AdsorptionReaction::getMolarFraction(x_react_, M_react_, M_carrier_);

    p_r_g_ = std::max(mol_frac_react * p_gas_, 1.0e-3); // avoid illdefined log
    setChemicalEquilibrium();
    const double dXdt = CaHydration();
    qR_ = (rho_up - rho_low) * dXdt;
}

// determine equilibrium temperature and pressure according to van't Hoff
void ReactionCaOH2::setChemicalEquilibrium()
{
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;

    X_D_ = (rho_s_ - rho_up - tol_rho_)/(rho_low - rho_up - 2.0*tol_rho_) ;
    X_D_ = (X_D_ < 0.5) ? std::max(tol_l_,X_D_) : std::min(X_D_,tol_u_); // constrain to interval [tol_l;tol_u]

    X_H_ = 1.0 - X_D_;

    // calculate equilibrium
    // using the p_eq to calculate the T_eq - Clausius-Clapeyron
    T_eq_ = (reaction_enthalpy_/R) / ((reaction_entropy_/R) + std::log(p_r_g_)); // unit of p in bar
    // Alternative: Use T_s as T_eq and calculate p_eq - for Schaube kinetics
    p_eq_ = std::exp((reaction_enthalpy_/R)/T_s_ - (reaction_entropy_/R));
}


double ReactionCaOH2::CaHydration()
{
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
    double dXdt;
        // step 3, calculate dX/dt
#ifdef SIMPLE_KINETICS
    if ( T_s < T_eq ) // hydration - simple model
#else
    if ( p_r_g_ > p_eq_ ) // hydration - Schaube model
#endif
    {
        //X_H = max(tol_l,X_H); //lower tolerance to avoid oscillations at onset of hydration reaction. Set here so that no residual reaction rate occurs at end of hydration.
#ifdef SIMPLE_KINETICS // this is from P. Schmidt
        dXdt = -1.0*(1.0-X_H) * (T_s - T_eq) / T_eq * 0.2 * conversion_rate::x_react;
#else //this is from Schaube
        if (X_H_ == tol_u_ || rho_s_ == rho_up)
        {
            dXdt = 0.0;
        }
        else if ((T_eq_ - T_s_) >= 50.0)
        {
            dXdt = 13945.0 * exp(-89486.0/R/T_s_) * std::pow(p_r_g_/p_eq_ - 1.0,0.83) * 3.0 * (X_D_) * std::pow(-1.0*log(X_D_),0.666);
        }
        else
        {
            dXdt = 1.0004e-34 * exp(5.3332e4 / T_s_) * std::pow(p_r_g_, 6.0) *
                   (X_D_);
        }
#endif
    }
    else // dehydration
    {
        //X_D = max(tol_l,X_D); //lower tolerance to avoid oscillations at onset of dehydration reaction. Set here so that no residual reaction rate occurs at end of dehydration.
#ifdef SIMPLE_KINETICS // this is from P. Schmidt
        dXdt = -1.0* (1.0-X_D) * (T_s - T_eq) / T_eq * 0.05;
#else
        if (X_D_ == tol_u_ || rho_s_ == rho_low)
        {
            dXdt = 0.0;
        }
        else if (X_D_ < 0.2)
        {
            dXdt = -1.9425e12 * exp( -1.8788e5/R/T_s_ ) * std::pow(1.0-p_r_g_/p_eq_,3.0)*(X_H_);
        }
        else
        {
            dXdt = -8.9588e9 * exp(-1.6262e5 / R / T_s_) *
                   std::pow(1.0 - p_r_g_ / p_eq_, 3.0) * 2.0 *
                   std::pow(X_H_, 0.5);
        }
#endif
    }
    return dXdt;
}

}  // namespace Adsorption
