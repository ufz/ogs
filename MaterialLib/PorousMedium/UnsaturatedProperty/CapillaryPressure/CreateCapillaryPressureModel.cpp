/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "BrooksCoreyCapillaryPressureSaturation.h"
#include "CapillaryPressureSaturation.h"
#include "CapillaryPressureSaturationCurve.h"
#include "VanGenuchtenCapillaryPressureSaturation.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
    \param config ConfigTree object which contains the input data
                  including `<type>BrooksCorey</type>`
                  and it has a tag of `<capillary_pressure>`
*/
static std::unique_ptr<CapillaryPressureSaturation> createBrooksCorey(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__capillary_pressure__type}
    config.checkConfigParameter("type", "BrooksCorey");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrooksCorey__pd}
    const auto pd = config.getConfigParameter<double>("pd");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrooksCorey__sr}
    const auto Sr = config.getConfigParameter<double>("sr");

    double Sg_r = 0.0;
    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrooksCorey__sg_r}
    if (auto const Sg_r_ptr = config.getConfigParameterOptional<double>("sg_r"))
    {
        DBUG(
            "Using value %g for nonwetting phase residual saturation in "
            "capillary pressure model.",
            (*Sg_r_ptr));
        Sg_r = *Sg_r_ptr;
    }
    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrooksCorey__smax}
    const auto Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrooksCorey__m}
    const auto m = config.getConfigParameter<double>("m");
    if (m < 1.0)  // m >= 1
    {
        OGS_FATAL(
            "The exponent parameter of BrooksCorey capillary pressure "
            "saturation model, m, must not be smaller than 1");
    }
    //! \ogs_file_param{material__porous_medium__capillary_pressure__BrooksCorey__pc_max}
    const auto Pc_max = config.getConfigParameter<double>("pc_max");

    return std::make_unique<BrooksCoreyCapillaryPressureSaturation>(
        pd, Sr, Sg_r, Smax, m, Pc_max);
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
    const auto pd = config.getConfigParameter<double>("pd");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__sr}
    const auto Sr = config.getConfigParameter<double>("sr");

    double Sg_r = 0.0;
    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__sg_r}
    if (auto const Sg_r_ptr = config.getConfigParameterOptional<double>("sg_r"))
    {
        DBUG(
            "Using value %g for nonwetting phase residual saturation in "
            "capillary pressure model.",
            (*Sg_r_ptr));
        Sg_r = *Sg_r_ptr;
    }

    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__smax}
    const auto Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__m}
    const auto m = config.getConfigParameter<double>("m");
    if (m < 0. || m > 1.0)
    {
        OGS_FATAL(
            "The exponent parameter of van Genuchten capillary pressure "
            "saturation model, m, must be in an interval of [0, 1]");
    }
    //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__pc_max}
    const auto Pc_max = config.getConfigParameter<double>("pc_max");

    bool has_regularized = false;
    if (auto const has_regularized_conf =
            //! \ogs_file_param{material__porous_medium__capillary_pressure__vanGenuchten__has_regularized}
        config.getConfigParameterOptional<bool>("has_regularized"))
    {
        DBUG("capillary pressure model: %s",
             (*has_regularized_conf) ? "true" : "false");
        has_regularized = *has_regularized_conf;
    }
    return std::make_unique<VanGenuchtenCapillaryPressureSaturation>(
        pd, Sr, Sg_r, Smax, m, Pc_max, has_regularized);
}

std::unique_ptr<CapillaryPressureSaturation> createCapillaryPressureModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__capillary_pressure__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "BrooksCorey")
    {
        return createBrooksCorey(config);
    }
    if (type == "vanGenuchten")
    {
        return createVanGenuchten(config);
    }
    if (type == "Curve")
    {
        //! \ogs_file_param{material__porous_medium__capillary_pressure__type}
        config.checkConfigParameter("type", "Curve");

        //! \ogs_file_param{material__porous_medium__capillary_pressure__Curve__curve}
        auto const& curve_config = config.getConfigSubtree("curve");

        auto curve = MathLib::createPiecewiseLinearCurve<
            MathLib::PiecewiseLinearMonotonicCurve>(curve_config);

        return std::make_unique<CapillaryPressureSaturationCurve>(
            std::move(curve));
    }

    OGS_FATAL(
        "The capillary pressure saturation models %s are unavailable.\n"
        "The available types are: \n\tBrooksCorey, \n\tvanGenuchten,",
        "\n\tCurve.\n",
        type.data());
}

}  // end namespace
}  // end namespace
