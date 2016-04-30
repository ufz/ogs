#pragma once

#include "Adsorption.h"

namespace Adsorption
{

class DensityConst : public AdsorptionReaction
{
public:
    double get_adsorbate_density(const double T_Ads) const;
    double get_alphaT(const double T_Ads) const;
    double characteristic_curve(const double A) const;
    double d_characteristic_curve(const double A) const;
};

}
