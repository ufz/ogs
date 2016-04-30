#ifndef ADSORPTION_H
#define ADSORPTION_H

#include <array>
#include <cmath>

#include "Reaction.h"

namespace Adsorption
{

const double GAS_CONST = 8.3144621;

const double M_N2  = 0.028;
const double M_H2O = 0.018;

class Adsorption : public Reaction
{
public:
    // TODO [CL] move those three methods to water properties class
    static double get_evaporation_enthalpy(const double T_Ads);
    static double get_equilibrium_vapour_pressure(const double T_Ads);
    static double get_specific_heat_capacity(const double T_Ads); // TODO [CL] why unused?

    static double get_molar_fraction(double xm, double M_this, double M_other);
    static double get_mass_fraction(double xn, double M_this, double M_other);
    static double d_molar_fraction(double xm, double M_this, double M_other);

    static double get_loading(const double rho_curr, const double rho_dry);

// non-virtual members
    double get_equilibrium_loading(const double p_Ads, const double T_Ads, const double M_Ads)
    const override;

// virtual members:
    virtual ~Adsorption() = default;

    virtual double get_enthalpy(const double p_Ads, const double T_Ads, const double M_Ads) const override;
    virtual double get_reaction_rate(const double p_Ads, const double T_Ads,
                                     const double M_Ads, const double loading) const override;
    /**
     * @brief get_d_reaction_rate
     * @param p_Ads
     * @param T_ads
     * @param M_Ads
     * @param loading
     * @param dqdr array containing the differentials wrt: p, T, C
     */
    virtual void get_d_reaction_rate(const double p_Ads, const double T_Ads,
                                     const double M_Ads, const double loading,
                                     std::array<double, 3>& dqdr) const;

protected:
    virtual double get_adsorbate_density(const double T_Ads) const = 0;
    virtual double get_alphaT(const double T_Ads) const = 0;
    virtual double characteristic_curve(const double A) const = 0;
    virtual double d_characteristic_curve(const double A) const = 0;

private:
// non-virtual members
    double get_potential(const double p_Ads, const double T_Ads, const double M_Ads) const;
    double get_entropy(const double T_Ads, const double A) const;
};


inline double curve_polyfrac(const double* coeffs, const double x)
{
    return ( coeffs[0] + coeffs[2] * x + coeffs[4] * pow(x,2) + coeffs[6] * pow(x,3) )
            / ( 1.0 + coeffs[1] * x + coeffs[3] * pow(x,2) + coeffs[5] * pow(x,3) );

    // Apparently, even pow and std::pow are different
    // return ( coeffs[0] + coeffs[2] * x + coeffs[4] * std::pow(x,2) + coeffs[6] * std::pow(x,3) )
    //         / ( 1.0 + coeffs[1] * x + coeffs[3] * std::pow(x,2) + coeffs[5] * std::pow(x,3) );

    // Analytically the same, but numerically quite different
    // return ( coeffs[0] + x * ( coeffs[2] + x * (coeffs[4] + x * coeffs[6] ) ) )
    //         / ( 1.0 + x * ( coeffs[1] + x * (coeffs[3] + x * coeffs[5] ) ) );

    // Analytically the same, but numerically quite different
    // return ( coeffs[0] + x * coeffs[2] + x*x * coeffs[4] + x*x*x * coeffs[6] )
    //        / ( 1.0 + x * coeffs[1] + x*x * coeffs[3] + x*x*x * coeffs[5] );
}

inline double d_curve_polyfrac(const double* coeffs, const double x)
{
    const double x2 = x*x;
    const double x3 = x2*x;
    const double u  = coeffs[0] + coeffs[2] * x +     coeffs[4] * x2 +     coeffs[6] * x3;
    const double du =             coeffs[2]     + 2.0*coeffs[4] * x  + 3.0*coeffs[6] * x2;
    const double v  = 1.0 + coeffs[1] * x +     coeffs[3] * x2 +     coeffs[5] * x3;
    const double dv =       coeffs[1]     + 2.0*coeffs[3] * x  + 3.0*coeffs[5] * x2;

    return (du*v - u*dv) / v / v;
}

}

#endif // ADSORPTION_H
