/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreateRelativePermeabilityModel.cpp
 *
 * Created on November 2, 2016, 11:43 AM
 */

#include "CreateRelativePermeabilityModel.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "MathLib/Curve/CreatePiecewiseLinearCurve.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "NonWettingPhaseBrooksCoreyOilGas.h"
#include "NonWettingPhaseVanGenuchten.h"
#include "RelativePermeability.h"
#include "RelativePermeabilityCurve.h"
#include "WettingPhaseBrooksCoreyOilGas.h"
#include "WettingPhaseVanGenuchten.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
    \param config ConfigTree object which contains the input data
                  including `<type>WettingPhaseVanGenuchten</type>`
                  and it has a tag of `<relative_permeability>`
*/
std::unique_ptr<RelativePermeability> createWettingPhaseVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__relative_permeability__type}
    config.checkConfigParameter("type", "WettingPhaseVanGenuchten");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseVanGenuchten__sr}
    const auto Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseVanGenuchten__smax}
    const auto Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseVanGenuchten__m}
    const auto m = config.getConfigParameter<double>("m");
    if (m < 0. || m > 1.0)
    {
        OGS_FATAL(
            "The exponent parameter of WettingPhaseVanGenuchten relative\n"
            " permeability model, m, must be in an interval of [0, 1]");
    }
    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseVanGenuchten__krel_min}
    const auto krel_min = config.getConfigParameter<double>("krel_min");

    return std::make_unique<WettingPhaseVanGenuchten>(Sr, Smax, m, krel_min);
}

/**
    \param config ConfigTree object which contains the input data
                  including `<type>NonWettingPhaseVanGenuchten</type>`
                  and it has a tag of `<relative_permeability>`
*/
std::unique_ptr<RelativePermeability> createNonWettingPhaseVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__relative_permeability__type}
    config.checkConfigParameter("type", "NonWettingPhaseVanGenuchten");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__sr}
    const auto Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__smax}
    const auto Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__m}
    const auto m = config.getConfigParameter<double>("m");
    if (m < 0. || m > 1.0)
    {
        OGS_FATAL(
            "The exponent parameter of NonWettingPhaseVanGenuchten relative\n"
            " permeability model, m, must be in an interval of [0, 1]");
    }

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__krel_min}
    const auto krel_min = config.getConfigParameter<double>("krel_min");

    return std::make_unique<NonWettingPhaseVanGenuchten>(Sr, Smax, m, krel_min);
}

/**
    \param config ConfigTree object which contains the input data
                  including `<type>WettingPhaseBrooksCoreyOilGas</type>`
                  and it has a tag of `<relative_permeability>`
*/
std::unique_ptr<RelativePermeability> createWettingPhaseBrooksCoreyOilGas(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__relative_permeability__type}
    config.checkConfigParameter("type", "WettingPhaseBrooksCoreyOilGas");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseBrooksCoreyOilGas__sr}
    const auto Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseBrooksCoreyOilGas__smax}
    const auto Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseBrooksCoreyOilGas__m}
    const auto m = config.getConfigParameter<double>("m");
    if (m < 1.0)  // m >= 1
    {
        OGS_FATAL(
            "The exponent parameter of WettingPhaseBrooksCoreyOilGas\n"
            "relative permeability model, m, must not be smaller than 1");
    }

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseBrooksCoreyOilGas__krel_min}
    const auto krel_min = config.getConfigParameter<double>("krel_min");

    return std::make_unique<WettingPhaseBrooksCoreyOilGas>(
        Sr, Smax, m, krel_min);
}

/**
    \param config ConfigTree object which contains the input data
                  including `<type>NonWettingPhaseBrooksCoreyOilGas</type>`
                  and it has a tag of `<relative_permeability>`
*/
std::unique_ptr<RelativePermeability> createNonWettingPhaseBrooksCoreyOilGas(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__relative_permeability__type}
    config.checkConfigParameter("type", "NonWettingPhaseBrooksCoreyOilGas");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseBrooksCoreyOilGas__sr}
    const auto Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseBrooksCoreyOilGas__smax}
    const auto Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseBrooksCoreyOilGas__m}
    const auto m = config.getConfigParameter<double>("m");
    if (m < 1.0)  // m >= 1
    {
        OGS_FATAL(
            "The exponent parameter of NonWettingPhaseBrooksCoreyOilGas\n"
            "relative permeability model, m, must not be smaller than 1");
    }

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseBrooksCoreyOilGas__krel_min}
    const auto krel_min = config.getConfigParameter<double>("krel_min");

    return std::make_unique<NonWettingPhaseBrooksCoreyOilGas>(
        Sr, Smax, m, krel_min);
}

std::unique_ptr<RelativePermeability> createRelativePermeabilityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__relative_permeability__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "WettingPhaseVanGenuchten")
    {
        return createWettingPhaseVanGenuchten(config);
    }
    if (type == "NonWettingPhaseVanGenuchten")
    {
        return createNonWettingPhaseVanGenuchten(config);
    }
    if (type == "WettingPhaseBrooksCoreyOilGas")
    {
        return createWettingPhaseBrooksCoreyOilGas(config);
    }
    if (type == "NonWettingPhaseBrooksCoreyOilGas")
    {
        return createNonWettingPhaseBrooksCoreyOilGas(config);
    }
    else if (type == "Curve")
    {
        //! \ogs_file_param{material__porous_medium__relative_permeability__type}
        config.checkConfigParameter("type", "Curve");

        //! \ogs_file_param{material__porous_medium__relative_permeability__Curve__curve}
        auto const& curve_config = config.getConfigSubtree("curve");

        auto curve = MathLib::createPiecewiseLinearCurve<MathLib
                              ::PiecewiseLinearInterpolation>(curve_config);
        return std::make_unique<RelativePermeabilityCurve>(std::move(curve));
    }
    else
    {
        OGS_FATAL(
            "The relative permeability model %s is unavailable.\n"
            "The available models are:"
            "\n\tWettingPhaseVanGenuchten,"
            "\n\tNonWettingPhaseVanGenuchten,"
            "\n\tWettingPhaseBrooksCoreyOilGas,"
            "\n\tNonWettingPhaseBrooksCoreyOilGas,",
            "\n\tCurve.\n",
            type.data());
    }
}

}  // end namespace
}  // end namespace
