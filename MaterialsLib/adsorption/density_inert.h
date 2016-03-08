#pragma once

#include "adsorption.h"

namespace Ads
{

class DensityInert final : public Adsorption
{
public:
	double get_enthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
						const double /*M_Ads*/) const override
	{
		return 0.0;
	}
	double get_reaction_rate(const double /*p_Ads*/, const double /*T_Ads*/,
							 const double /*M_Ads*/, const double /*loading*/
							 ) const override
	{
		return 0.0;
	}

protected:
	double get_adsorbate_density(const double /*T_Ads*/) const override
	{
		return 1.0;
	}
	double get_alphaT(const double /*T_Ads*/) const override
	{
		return 0.0;
	}
	double characteristic_curve(const double A) const override
	{
		return A;
	}
	double d_characteristic_curve(const double /*A*/) const override
	{
		return 1.0;
	}
};

}
