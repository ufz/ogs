/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Adsorption.h"

namespace Adsorption
{

class DensityMette : public AdsorptionReaction
{
public:
    double getAdsorbateDensity(const double T_Ads) const;
    double getAlphaT(const double T_Ads) const;
    double characteristicCurve(const double A) const;
    double dCharacteristicCurve(const double A) const;
};

}
