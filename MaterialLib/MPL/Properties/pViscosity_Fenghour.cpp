/**
 * \author Norbert Grunwald
 * \date   Sep 21, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pViscosity_Fenghour.h"
#include "../mpComponent.h"

namespace MaterialPropertyLib
{
ViscosityCO2Fenghour::ViscosityCO2Fenghour(Medium*) : _component(0)
{
    notImplemented("ViscosityCO2Fenghour", "Medium");
};

ViscosityCO2Fenghour::ViscosityCO2Fenghour(Phase*) : _component(0)
{
    notImplemented("ViscosityCO2Fenghour", "Phase");
};

ViscosityCO2Fenghour::ViscosityCO2Fenghour(Component* c)
    : _component(c),
      _a{0.235156, -0.491266, 0.05211155, 0.05347906, -0.01537102},
      _d{0.4071119e-02, 0.7198037e-04, 0.2411697e-16, 0.2971072e-22,
         -0.1627888e-22} {};

PropertyDataType ViscosityCO2Fenghour::value(VariableArray const& vars)
{
    if (isUpdated())
        return _value;

    const double temperature = getScalar(vars[T]);
    const double epsilon_k = 251.196;  // I think it might be ok to hard-code
    // such rare material constants; It is not used anywhere else, and this
    // method here is specifically designed for carbon doxide...
    const double T_red = temperature / epsilon_k;
    const double rho = getScalar(_component->property(density), vars);

    // some powers of rho, just in order to avoid expensive pow methods
    const double rho_pow_2 = rho * rho;
    const double rho_pow_6 = rho_pow_2 * rho_pow_2 * rho_pow_2;
    const double rho_pow_8 = rho_pow_6 * rho_pow_2;

    double psi(0.);

    for (std::size_t i = 0; i < _a.size(); i++)
        psi += _a[i] * std::pow(std::log(T_red), i);

    const double eta_0 = 1.00697 * std::sqrt(T_red) / std::exp(psi);
    const double eta_d = _d[0] * rho + _d[1] * rho_pow_2 +
                         _d[2] * rho_pow_6 / T_red + _d[3] * rho_pow_8 +
                         _d[4] * rho_pow_8 / T_red;
    const double eta = (eta_0 + eta_d) / 1.e6; // conversion MPa*s in Pa*s

    _value = eta;
    return eta;
}

}  // MaterialPropertyLib
