/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionCaOH2.h"
#include <cassert>
#include <logog/include/logog.hpp>
#include "BaseLib/Error.h"
#include "MaterialLib/PhysicalConstant.h"
#include "Adsorption.h"

namespace Adsorption
{

const double ReactionCaOH2::_reaction_enthalpy = -1.12e+05;
const double ReactionCaOH2::_reaction_entropy  = -143.5;
const double ReactionCaOH2::_M_carrier =
    MaterialLib::PhysicalConstant::MolarMass::N2;
const double ReactionCaOH2::_M_react =
    MaterialLib::PhysicalConstant::MolarMass::Water;

const double ReactionCaOH2::_tol_l   = 1e-4;
const double ReactionCaOH2::_tol_u   = 1.0 - 1e-4;
const double ReactionCaOH2::_tol_rho = 0.1;

const double ReactionCaOH2::rho_low = 1656.0;
const double ReactionCaOH2::rho_up = 2200.0;


double
ReactionCaOH2::getEnthalpy(const double, const double, const double) const
{
    return - _reaction_enthalpy/_M_react;
}

double
ReactionCaOH2::getReactionRate(const double, const double, const double, const double) const
{
    OGS_FATAL("get_reaction_rate do not call directly");
}


double ReactionCaOH2::getReactionRate(double const solid_density)
{
    _rho_s = solid_density;
    calculateQR();
    return _qR;
}

void ReactionCaOH2::updateParam(
        double T_solid,
        double p_gas,
        double x_react,
        double rho_s_initial)
{
    _T_s     = T_solid;
    _p_gas = p_gas / 1e5; // convert Pa to bar
    _x_react = x_react;
    _rho_s   = rho_s_initial;
}

void ReactionCaOH2::calculateQR()
{
    // Convert mass fraction into mole fraction
    const double mol_frac_react = AdsorptionReaction::getMolarFraction(_x_react, _M_react, _M_carrier);

    _p_r_g = std::max(mol_frac_react * _p_gas, 1.0e-3); // avoid illdefined log
    setChemicalEquilibrium();
    const double dXdt = CaHydration();
    _qR = (rho_up - rho_low) * dXdt;
}

// determine equilibrium temperature and pressure according to van't Hoff
void ReactionCaOH2::setChemicalEquilibrium()
{
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;

    _X_D = (_rho_s - rho_up - _tol_rho)/(rho_low - rho_up - 2.0*_tol_rho) ;
    _X_D = (_X_D < 0.5) ? std::max(_tol_l,_X_D) : std::min(_X_D,_tol_u); // constrain to interval [tol_l;tol_u]

    _X_H = 1.0 - _X_D;

    // calculate equilibrium
    // using the p_eq to calculate the T_eq - Clausius-Clapeyron
    _T_eq = (_reaction_enthalpy/R) / ((_reaction_entropy/R) + std::log(_p_r_g)); // unit of p in bar
    // Alternative: Use T_s as T_eq and calculate p_eq - for Schaube kinetics
    _p_eq = std::exp((_reaction_enthalpy/R)/_T_s - (_reaction_entropy/R));
}


double ReactionCaOH2::CaHydration()
{
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
    double dXdt;
        // step 3, calculate dX/dt
#ifdef SIMPLE_KINETICS
    if ( T_s < T_eq ) // hydration - simple model
#else
    if ( _p_r_g > _p_eq ) // hydration - Schaube model
#endif
    {
        //X_H = max(tol_l,X_H); //lower tolerance to avoid oscillations at onset of hydration reaction. Set here so that no residual reaction rate occurs at end of hydration.
#ifdef SIMPLE_KINETICS // this is from P. Schmidt
        dXdt = -1.0*(1.0-X_H) * (T_s - T_eq) / T_eq * 0.2 * conversion_rate::x_react;
#else //this is from Schaube
        if (_X_H == _tol_u || _rho_s == rho_up)
            dXdt = 0.0;
        else if ( (_T_eq-_T_s) >= 50.0)
            dXdt = 13945.0 * exp(-89486.0/R/_T_s) * std::pow(_p_r_g/_p_eq - 1.0,0.83) * 3.0 * (_X_D) * std::pow(-1.0*log(_X_D),0.666);
        else
            dXdt = 1.0004e-34 * exp(5.3332e4/_T_s) * std::pow(_p_r_g, 6.0) * (_X_D);
#endif
    }
    else // dehydration
    {
        //X_D = max(tol_l,X_D); //lower tolerance to avoid oscillations at onset of dehydration reaction. Set here so that no residual reaction rate occurs at end of dehydration.
#ifdef SIMPLE_KINETICS // this is from P. Schmidt
        dXdt = -1.0* (1.0-X_D) * (T_s - T_eq) / T_eq * 0.05;
#else
        if (_X_D == _tol_u || _rho_s == rho_low)
            dXdt = 0.0;
        else if (_X_D < 0.2)
            dXdt = -1.9425e12 * exp( -1.8788e5/R/_T_s ) * std::pow(1.0-_p_r_g/_p_eq,3.0)*(_X_H);
        else
            dXdt = -8.9588e9 * exp( -1.6262e5/R/_T_s ) * std::pow(1.0-_p_r_g/_p_eq,3.0)*2.0 * std::pow(_X_H, 0.5);
#endif
    }
    return dXdt;
}

}
