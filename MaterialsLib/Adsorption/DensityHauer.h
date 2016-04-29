#pragma once

#include "Adsorption.h"
#include "DensityCook.h"

namespace Adsorption
{

class DensityHauer : public Adsorption
{
public:
	double get_adsorbate_density(const double T_Ads) const;
	double get_alphaT(const double T_Ads) const;
	double characteristic_curve(const double A) const;
	double d_characteristic_curve(const double A) const;
};

inline double rho_water_Hauer(const double T_Ads)
{
	// data like in python script
	const double T0 = 283.15, rho0 = rho_water_Dean(T0), alpha0 = 3.781e-4; //K; kg/m^3; 1/K

	return rho0 * (1. - alpha0 * (T_Ads-T0)); //in kg/m^3
}

}
