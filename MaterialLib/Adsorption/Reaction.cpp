/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <logog/include/logog.hpp>

#include "BaseLib/StringTools.h"

#include "Reaction.h"

#include "Density100MPa.h"
#include "DensityConst.h"
#include "DensityCook.h"
#include "DensityDubinin.h"
#include "DensityHauer.h"
#include "DensityLegacy.h"
#include "DensityMette.h"
#include "DensityNunez.h"

#include "ReactionCaOH2.h"
#include "ReactionInert.h"
#include "ReactionSinusoidal.h"


namespace Adsorption
{

std::unique_ptr<Reaction>
Reaction::
newInstance(BaseLib::ConfigTree const& conf)
{
    //! \ogs_file_param{material__adsorption__reaction__type}
    auto const type = conf.getConfigParameter<std::string>("type");

    if (type == "Z13XBF")
        return std::make_unique<DensityLegacy>();
    if (type == "Z13XBF_100MPa")
        return std::make_unique<Density100MPa>();
    if (type == "Z13XBF_Const")
        return std::make_unique<DensityConst>();
    if (type == "Z13XBF_Cook")
        return std::make_unique<DensityCook>();
    if (type == "Z13XBF_Dubinin")
        return std::make_unique<DensityDubinin>();
    if (type == "Z13XBF_Hauer")
        return std::make_unique<DensityHauer>();
    if (type == "Z13XBF_Mette")
        return std::make_unique<DensityMette>();
    if (type == "Z13XBF_Nunez")
        return std::make_unique<DensityNunez>();
    if (type == "Inert")
        return std::make_unique<ReactionInert>();
    if (type == "Sinusoidal")
        return std::make_unique<ReactionSinusoidal>(conf);
    if (type == "CaOH2")
        return std::make_unique<ReactionCaOH2>(conf);

    OGS_FATAL("Unknown reactive system: %s.", type.c_str());

    return nullptr;
}

} // namespace Adsorption
