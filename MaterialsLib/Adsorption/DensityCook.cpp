#include "DensityCook.h"

namespace
{

// NaX_Dean_polyfrac_CC.pickle
// date extracted 2015-06-23 15:38:35 file mtime 2015-06-23 15:19:42
const double c[] = {
	0.3632627555646154,		/* a0 */
	-0.0014090624975800715,	/* a1 */
	-0.0007717609035743321,	/* a2 */
	5.03634836561135e-09,	/* a3 */
	5.478509959282738e-07,	/* a4 */
	6.36458510620815e-10,	/* a5 */
	-1.037977321231462e-10	/* a6 */
};

}

namespace Ads
{

double DensityCook::get_adsorbate_density(const double T_Ads) const
{
	return rho_water_Dean(T_Ads);
}


//Thermal expansivity model for water found in the works of Hauer
double DensityCook::get_alphaT(const double T_Ads) const
{
	return alphaT_water_Dean(T_Ads);
}


//Characteristic curve. Return W (A)
double DensityCook::characteristic_curve(const double A) const
{
	double W = curve_polyfrac(c, A); //cm^3/g

	if (W < 0.0) {
		W = 0.0; // TODO [CL] debug output
	}

	return W/1.e3; //m^3/kg
}

double DensityCook::d_characteristic_curve(const double A) const
{
	return d_curve_polyfrac(c, A);
}

}
