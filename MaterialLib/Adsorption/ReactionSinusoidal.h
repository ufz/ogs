/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALSLIB_ADSORPTION_REACTIONSINUSOIDAL_H
#define MATERIALSLIB_ADSORPTION_REACTIONSINUSOIDAL_H

#include <logog/include/logog.hpp>

#include "Reaction.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/StringTools.h"

namespace Adsorption
{

class ReactionSinusoidal final : public Reaction
{
public:
    explicit ReactionSinusoidal(BaseLib::ConfigTree const& conf) :
        //! \ogs_file_param{material__adsorption__reaction__Sinusoidal__reaction_enthalpy}
        _enthalpy(conf.getConfigParameter<double>("reaction_enthalpy"))
    {
    }

    double getEnthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
                        const double /*M_Ads*/) const override
    {
        return _enthalpy;
    }

    double getReactionRate(const double /*p_Ads*/, const double /*T_Ads*/, const double /*M_Ads*/,
                             const double /*loading*/) const override
    {
        OGS_FATAL("Method getReactionRate() should never be called directly");
    }

private:
    double _enthalpy;
};

}
#endif // MATERIALSLIB_ADSORPTION_REACTIONSINUSOIDAL_H
