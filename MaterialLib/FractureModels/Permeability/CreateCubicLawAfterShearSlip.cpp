/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateCubicLawAfterShearSlip.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "CubicLawAfterShearSlip.h"

namespace MaterialLib::Fracture::Permeability
{
std::unique_ptr<Permeability> createCubicLawAfterShearSlip(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fracture_properties__permeability_model__type}
    config.checkConfigParameter("type", "CubicLawAfterShearSlip");

    auto const initial_creation_aperture =
        //! \ogs_file_param{material__fracture_properties__permeability_model__CubicLawAfterShearSlip__initial_creation_aperture}
        config.getConfigParameter<double>("initial_creation_aperture");

    auto const minimum_permeability =
        //! \ogs_file_param{material__fracture_properties__permeability_model__CubicLawAfterShearSlip__minimum_permeability}
        config.getConfigParameter<double>("minimum_permeability");

    auto const aperture_threshold =
        //! \ogs_file_param{material__fracture_properties__permeability_model__CubicLawAfterShearSlip__aperture_threshold}
        config.getConfigParameter<double>("aperture_threshold");

    if (minimum_permeability < 0)
    {
        OGS_FATAL(
            "The minimum_permeability parameter must be non-negative. Given "
            "value %g.",
            minimum_permeability);
    }

    return std::make_unique<CubicLawAfterShearSlip>(
        initial_creation_aperture, minimum_permeability, aperture_threshold);
}
}  // namespace MaterialLib::Fracture::Permeability
