/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreateCapillaryPressureModel.cpp
 *
 * Created on November 1, 2016, 10:06 AM
 */

#include "CreateCapillaryPressureModel.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "MathLib/Curve/CreatePiecewiseLinearCurve.h"
#include "MathLib/Curve/PiecewiseLinearMonotonicCurve.h"

#include "BrookCoreyCapillaryPressureSaturation.h"
#include "CapillaryPressureSaturation.h"
#include "CapillaryPressureSaturationCurve.h"
#include "VanGenuchtenCapillaryPressureSaturation.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
    \param config ConfigTree object which contains the input data
                  including `<type>BrookCorey</type>`
                  and it has a tag of `<capillary_pressure>`
*/
static std::unique_ptr<CapillaryPressureSaturation> createBrookCorey(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__capillary_pressure__type}
    config.checkConfigParameter("type", "BrookCorey");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrookCorey__pd}
    const double pd = config.getConfigParameter<double>("pd");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrookCorey__sr}
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrookCorey__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrookCorey__m}
    const double m = config.getConfigParameter<double>("m");
    if (m < 1.0)  // m >= 1
    {
        OGS_FATAL(
            "The exponent parameter of BrookCorey capillary pressure "
            "saturation model, m, must not be smaller than 1");
    }
    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrookCorey__pc_max}
    const double Pc_max = config.getConfigParameter<double>("pc_max");

    return std::unique_ptr<CapillaryPressureSaturation>(
        new BrookCoreyCapillaryPressureSaturation(pd, Sr, Smax, m, Pc_max));
}

/**
    \param config ConfigTree object which contains the input data
                  including `<type>vanGenuchten</type>`
                  and it has a tag of `<capillary_pressure>`
*/
static std::unique_ptr<CapillaryPressureSaturation> createVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__capillary_pressure__type}
    config.checkConfigParameter("type", "vanGenuchten");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__pd}
    const double pd = config.getConfigParameter<double>("pd");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__sr}
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__m}
    const double m = config.getConfigParameter<double>("m");
    if (m < 0. || m > 1.0)
    {
        OGS_FATAL(
            "The exponent parameter of van Genuchten capillary pressure "
            "saturation model, m, must be in an interval of [0, 1]");
    }
    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__pc_max}
    const double Pc_max = config.getConfigParameter<double>("pc_max");

    return std::unique_ptr<CapillaryPressureSaturation>(
        new VanGenuchtenCapillaryPressureSaturation(pd, Sr, Smax, m, Pc_max));
}

std::unique_ptr<CapillaryPressureSaturation> createCapillaryPressureModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__capillary_pressure__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "BrookCorey")
    {
        return createBrookCorey(config);
    }
    else if (type == "vanGenuchten")
    {
        return createVanGenuchten(config);
    }
    else if (type == "Curve")
    {
        //! \ogs_file_param{material__porous_medium__capillary_pressure__type}
        config.checkConfigParameter("type", "Curve");

        //! \ogs_file_param{material__porous_medium__capillary_pressure__Curve__curve}
        auto const& curve_config = config.getConfigSubtree("curve");

        auto curve = MathLib::createPiecewiseLinearCurve<
            MathLib::PiecewiseLinearMonotonicCurve>(curve_config);

        return std::unique_ptr<CapillaryPressureSaturation>(
            new CapillaryPressureSaturationCurve(std::move(curve)));
    }
    else
    {
        OGS_FATAL(
            "The capillary pressure saturation models %s are unavailable.\n"
            "The available types are: \n\tBrookCorey, \n\tvanGenuchten,",
            "\n\tCurve.\n",
            type.data());
    }
}

}  // end namespace
}  // end namespace
