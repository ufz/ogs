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

#include "RelativePermeability.h"
#include "RelativePermeabilityCurve.h"
#include "WettingPhaseVanGenuchten.h"
#include "NonWettingPhaseVanGenuchten.h"
#include "WettingPhaseBrooksCoreyOilGas.h"
#include "NonWettingPhaseBrooksCoreyOilGas.h"

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
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseVanGenuchten__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseVanGenuchten__m}
    const double m = config.getConfigParameter<double>("m");
    if (m < 0. || m > 1.0)
    {
        OGS_FATAL(
            "The exponent parameter of WettingPhaseVanGenuchten relative\n"
            " permeability model, m, must be in an interval of [0, 1]");
    }
    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseVanGenuchten__krel_min}
    const double krel_min = config.getConfigParameter<double>("krel_min");

    return std::unique_ptr<RelativePermeability>(
        new WettingPhaseVanGenuchten(Sr, Smax, m, krel_min));
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
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__m}
    const double m = config.getConfigParameter<double>("m");
    if (m < 0. || m > 1.0)
    {
        OGS_FATAL(
            "The exponent parameter of NonWettingPhaseVanGenuchten relative\n"
            " permeability model, m, must be in an interval of [0, 1]");
    }

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__krel_min}
    const double krel_min = config.getConfigParameter<double>("krel_min");

    return std::unique_ptr<RelativePermeability>(
        new NonWettingPhaseVanGenuchten(Sr, Smax, m, krel_min));
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
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseBrooksCoreyOilGas__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseBrooksCoreyOilGas__m}
    const double m = config.getConfigParameter<double>("m");
    if (m < 1.0)  // m >= 1
    {
        OGS_FATAL(
            "The exponent parameter of WettingPhaseBrooksCoreyOilGas\n"
            "relative permeability model, m, must not be smaller than 1");
    }

    //! \ogs_file_param{material__porous_medium__relative_permeability__WettingPhaseBrooksCoreyOilGas__krel_min}
    const double krel_min = config.getConfigParameter<double>("krel_min");

    return std::unique_ptr<RelativePermeability>(
        new WettingPhaseBrooksCoreyOilGas(Sr, Smax, m, krel_min));
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
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseBrooksCoreyOilGas__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseBrooksCoreyOilGas__m}
    const double m = config.getConfigParameter<double>("m");
    if (m < 1.0)  // m >= 1
    {
        OGS_FATAL(
            "The exponent parameter of NonWettingPhaseBrooksCoreyOilGas\n"
            "relative permeability model, m, must not be smaller than 1");
    }

    //! \ogs_file_param{material__porous_medium__relative_permeability__NonWettingPhaseBrooksCoreyOilGas__krel_min}
    const double krel_min = config.getConfigParameter<double>("krel_min");

    return std::unique_ptr<RelativePermeability>(
        new NonWettingPhaseBrooksCoreyOilGas(Sr, Smax, m, krel_min));
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
    else if (type == "NonWettingPhaseVanGenuchten")
    {
        return createNonWettingPhaseVanGenuchten(config);
    }
    else if (type == "WettingPhaseBrooksCoreyOilGas")
    {
        return createWettingPhaseBrooksCoreyOilGas(config);
    }
    else if (type == "NonWettingPhaseBrooksCoreyOilGas")
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
        return std::unique_ptr<RelativePermeability>(
            new RelativePermeabilityCurve(std::move(curve)));
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
