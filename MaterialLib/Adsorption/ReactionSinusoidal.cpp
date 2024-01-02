/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ReactionSinusoidal.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

namespace Adsorption
{

ReactionSinusoidal::ReactionSinusoidal(BaseLib::ConfigTree const& conf)
    :  //! \ogs_file_param{material__adsorption__reaction__Sinusoidal__reaction_enthalpy}
      _enthalpy(conf.getConfigParameter<double>("reaction_enthalpy"))
{
}

double ReactionSinusoidal::getReactionRate(const double /*p_Ads*/,
                                           const double /*T_Ads*/,
                                           const double /*M_Ads*/,
                                           const double /*loading*/) const
{
    OGS_FATAL("Method getReactionRate() should never be called directly");
}

}  // namespace Adsorption
