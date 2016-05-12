/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Adsorption.h"

namespace Adsorption
{

class DensityLegacy : public AdsorptionReaction
{
public:
    double get_adsorbate_density(const double T_Ads) const;
    double get_alphaT(const double T_Ads) const;
    double characteristic_curve(const double A) const;
    double d_characteristic_curve(const double A) const;
};

}
