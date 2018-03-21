/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FracturePermeability.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

namespace MaterialLib
{
namespace Fracture
{
std::unique_ptr<ConstantPermeability> createConstantPermeability(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "ConstantPermeability");

    auto const permeability =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties__permeability}
        config.getConfigParameter<double>("value");

    if (permeability < 0)
        OGS_FATAL(
            "The permeability parameter must be non-negative. Given value %g.",
            permeability);

    return std::make_unique<ConstantPermeability>(permeability);
}

std::unique_ptr<ConstantHydraulicAperture> createConstantHydraulicAperture(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "ConstantHydraulicAperture");

    return std::make_unique<ConstantHydraulicAperture>();
}

std::unique_ptr<CubicLaw> createCubicLaw(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "CubicLaw");

    return std::make_unique<CubicLaw>();
}

std::unique_ptr<CubicLawAfterShearSlip> createCubicLawAfterShearSlip(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "CubicLawAfterShearSlip");

    auto const initial_creation_aperture =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties__initial_creation_aperture}
        config.getConfigParameter<double>("initial_creation_aperture");

    auto const minimum_permeability =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS_WITH_LIE__fracture_properties__minimum_permeability}
        config.getConfigParameter<double>("minimum_permeability");

    if (minimum_permeability < 0)
        OGS_FATAL(
            "The minimum_permeability parameter must be non-negative. Given "
            "value %g.",
            minimum_permeability);

    return std::make_unique<CubicLawAfterShearSlip>(initial_creation_aperture,
                                                    minimum_permeability);
}
}  // namespace Fracture
}  // namespace MaterialLib
