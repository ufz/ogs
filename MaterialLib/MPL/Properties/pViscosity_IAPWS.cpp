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

ViscosityWaterIAPWS::ViscosityWaterIAPWS(Component* c)
    : _component(c),
      _H{{1.67752, 2.20462, 0.6366564, -0.241605}},
      _h{{{0.5200940, 0.2225310, -0.2813780, 0.1619130, -0.0325372, 0., 0.},
          {0.0850895, 0.9991150, -0.9068510, 0.2573990, 0., 0., 0.},
          {-1.0837400, 1.8879700, -0.772479, 0., 0., 0., 0.},
          {-0.289555, 1.266130, -0.489837, 0., 0.0698452, 0., -0.00435673},
          {0., 0., -0.257040, 0., 0., 0.00872102, 0.},
          {0., 0.120573, 0., 0., 0., 0., -0.000593264}}} {};

PropertyDataType ViscosityWaterIAPWS::value(VariableArray const& vars)
{
    if (isUpdated())
        return _value;

    const double T_crit = getScalar(_component->property(critical_temperature));
    const double T_red = getScalar(vars[T]) / T_crit;
    const double rho_crit = getScalar(_component->property(critical_density));
    const double rho = getScalar(_component->property(density), vars);
    const double rho_red = rho / rho_crit;

    std::array<double, 3> sum{{0., 0., 0.}};

    for (unsigned i = 0; i < 4; ++i)
        sum[0] += _H[i] / std::pow(T_red, i);

    for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 7; ++j)
        {
            sum[2] += _h[i][j] * std::pow(rho_red - 1., j);
            sum[1] += (std::pow(1. / T_red - 1., i) * sum[2]);
            sum[2] = 0;
        }

    const double eta =
        100. * std::sqrt(T_red) / sum[0] * std::exp(rho_red * sum[1]) / 1.e6;

    _value = eta;
    return eta;
}

}  // MaterialPropertyLib
