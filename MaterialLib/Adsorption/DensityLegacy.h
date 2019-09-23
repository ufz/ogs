/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Adsorption.h"

namespace Adsorption
{

class DensityLegacy : public AdsorptionReaction
{
public:
    double getAdsorbateDensity(const double T_Ads) const override;
    double getAlphaT(const double T_Ads) const override;
    double characteristicCurve(const double A) const override;
    double dCharacteristicCurve(const double A) const override;
};

}  // namespace Adsorption
