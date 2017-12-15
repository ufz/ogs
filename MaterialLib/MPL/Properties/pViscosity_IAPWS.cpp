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

#include "pViscosity_IAPWS.h"
#include <cmath>
#include "../mpComponent.h"

namespace MaterialPropertyLib
{
ViscosityWaterIAPWS::ViscosityWaterIAPWS(Medium*) : _component(0)
{
    notImplemented("Viscosity_IAPWS", "Medium");
}

ViscosityWaterIAPWS::ViscosityWaterIAPWS(Phase*) : _component(0)
{
    notImplemented("Viscosity_IAPWS", "Phase");
}

ViscosityWaterIAPWS::ViscosityWaterIAPWS(Component* c) : _component(c)
{
}

static constexpr std::array<double, 4> H = {
    {1.67752, 2.20462, 0.6366564, -0.241605}};

static constexpr std::array<std::array<double, 7>, 6> h = {
    {{0.5200940, 0.2225310, -0.2813780, 0.1619130, -0.0325372, 0., 0.},
     {0.0850895, 0.9991150, -0.9068510, 0.2573990, 0., 0., 0.},
     {-1.0837400, 1.8879700, -0.772479, 0., 0., 0., 0.},
     {-0.289555, 1.266130, -0.489837, 0., 0.0698452, 0., -0.00435673},
     {0., 0., -0.257040, 0., 0., 0.00872102, 0.},
     {0., 0.120573, 0., 0., 0., 0., -0.000593264}}};

PropertyDataType ViscosityWaterIAPWS::value(VariableArray const& vars)
{
    if (isUpdated())
        return _value;

    const double T_crit = getScalar(_component->property(critical_temperature));
    const double T_red = getScalar(vars[T]) / T_crit;
    const double rho_crit = getScalar(_component->property(critical_density));
    const double rho = getScalar(_component->property(density), vars);
    const double rho_red = rho / rho_crit;

    double eta_0 = H.back();

    for (int i = H.size()-2; i >= 0; --i)
        eta_0 = eta_0/T_red + H[i];

    eta_0 = 100.*std::sqrt(T_red)/eta_0;

    double A = 1/T_red - 1.;
    double B = rho_red - 1.;

    double eta_1 = 0;

    const auto hi_max = h.size()-1;
    const auto hj_max = h[0].size()-1;

    for (int j=hj_max; j>=0; j--)
    {
        double a = h[hi_max][j];
        for (int i=hi_max-1; i>=0; i--)
            a = h[i][j] + A*a;
        eta_1 = a + B*eta_1;
    }

    eta_1 = std::exp(rho_red*eta_1);

    const double eta = eta_0 * eta_1 / 1.e6;

    _value = eta;
    return eta;
}

}  // MaterialPropertyLib
