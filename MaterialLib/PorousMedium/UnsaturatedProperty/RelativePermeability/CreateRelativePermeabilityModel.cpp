/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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
#include <sstream>
#include <string>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "RelativePermeability.h"
#include "ReletivePermeabilityCurve.h"
#include "WettingPhaseVanGenuchten.h"
#include "NonWettingPhaseVanGenuchten.h"
#include "WettingPhaseBrookCoreyOilGas.h"
#include "NonWettingPhaseBrookCoreyOilGas.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
    \param config ConfigTree object which contains the input data
                  including <type>WettingPhaseVanGenuchten</type>
                  and it has a tag of <relative_permeability>
*/
std::unique_ptr<RelativePermeability> createWettingPhaseVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type}
    config.checkConfigParameter("type", "WettingPhaseVanGenuchten");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__WettingPhaseVanGenuchten__sr}
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__WettingPhaseVanGenuchten__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__WettingPhaseVanGenuchten__m}
    const double m = config.getConfigParameter<double>("m");
    if (m > 1.0)  // m <= 1
    {
        OGS_FATAL(
            "The exponent parameter of WettingPhaseVanGenuchten relative\n"
            " permeability model, m, must be in an interval of [0, 1]");
    }
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__WettingPhaseVanGenuchten__m}
    const double krel_max = config.getConfigParameter<double>("krel_min");

    return std::unique_ptr<RelativePermeability>(
        new WettingPhaseVanGenuchten(Sr, Smax, m, krel_max));
}

/**
    \param config ConfigTree object which contains the input data
                  including <type>NonWettingPhaseVanGenuchten</type>
                  and it has a tag of <relative_permeability>
*/
std::unique_ptr<RelativePermeability> createNonWettingPhaseVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type}
    config.checkConfigParameter("type", "NonWettingPhaseVanGenuchten");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__sr}
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__m}
    const double m = config.getConfigParameter<double>("m");
    if (m > 1.0)  // m <= 1
    {
        OGS_FATAL(
            "The exponent parameter of NonWettingPhaseVanGenuchten relative\n"
            " permeability model, m, must be in an interval of [0, 1]");
    }
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__NonWettingPhaseVanGenuchten__m}
    const double krel_max = config.getConfigParameter<double>("krel_min");

    return std::unique_ptr<RelativePermeability>(
        new NonWettingPhaseVanGenuchten(Sr, Smax, m, krel_max));
}

/**
    \param config ConfigTree object which contains the input data
                  including <type>WettingPhaseBrookCoreyOilGas</type>
                  and it has a tag of <relative_permeability>
*/
std::unique_ptr<RelativePermeability> createWettingPhaseBrookCoreyOilGas(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type}
    config.checkConfigParameter("type", "WettingPhaseBrookCoreyOilGas");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__WettingPhaseBrookCoreyOilGas__sr}
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__WettingPhaseBrookCoreyOilGas__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__WettingPhaseBrookCoreyOilGas__m}
    const double m = config.getConfigParameter<double>("m");
    if (m < 1.0)  // m >= 1
    {
        OGS_FATAL(
            "The exponent parameter of WettingPhaseBrookCoreyOilGas\n"
            "relative permeability model, m, must not be smaller than 1");
    }
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__WettingPhaseBrookCoreyOilGas__m}
    const double krel_max = config.getConfigParameter<double>("krel_min");

    return std::unique_ptr<RelativePermeability>(
        new WettingPhaseBrookCoreyOilGas(Sr, Smax, m, krel_max));
}

/**
    \param config ConfigTree object which contains the input data
                  including <type>NonWettingPhaseBrookCoreyOilGas</type>
                  and it has a tag of <relative_permeability>
*/
std::unique_ptr<RelativePermeability> createNonWettingPhaseBrookCoreyOilGas(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type}
    config.checkConfigParameter("type", "NonWettingPhaseBrookCoreyOilGas");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__NonWettingPhaseBrookCoreyOilGas__sr}
    const double Sr = config.getConfigParameter<double>("sr");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__NonWettingPhaseBrookCoreyOilGas__smax}
    const double Smax = config.getConfigParameter<double>("smax");

    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__NonWettingPhaseBrookCoreyOilGas__m}
    const double m = config.getConfigParameter<double>("m");
    if (m < 1.0)  // m >= 1
    {
        OGS_FATAL(
            "The exponent parameter of NonWettingPhaseBrookCoreyOilGas\n"
            "relative permeability model, m, must not be smaller than 1");
    }
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__NonWettingPhaseBrookCoreyOilGas__m}
    const double krel_max = config.getConfigParameter<double>("krel_min");

    return std::unique_ptr<RelativePermeability>(
        new NonWettingPhaseBrookCoreyOilGas(Sr, Smax, m, krel_max));
}

std::unique_ptr<RelativePermeability> createRelativePermeabilityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "WettingPhaseVanGenuchten")
    {
        return createWettingPhaseVanGenuchten(config);
    }
    else if (type == "NonWettingPhaseVanGenuchten")
    {
        return createNonWettingPhaseVanGenuchten(config);
    }
    else if (type == "WettingPhaseBrookCoreyOilGas")
    {
        return createWettingPhaseBrookCoreyOilGas(config);
    }
    else if (type == "NonWettingPhaseBrookCoreyOilGas")
    {
        return createNonWettingPhaseBrookCoreyOilGas(config);
    }
    else if (type == "Curve")
    {
        config.checkConfigParameter("type", "Curve");

        std::vector<double> variables, values;
        //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__Curve__curve}
        auto const& curve_config = config.getConfigSubtree("curve");
        for (
            auto const& data_string :
            //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__Curve_curve__data}
            curve_config.getConfigParameterList<std::string>("data"))
        {
            std::stringstream ss(data_string);
            double var, val;
            ss >> var >> val;
            ss.clear();
            variables.push_back(var);
            values.push_back(val);
        }
        auto curve = std::unique_ptr<MathLib::PiecewiseLinearInterpolation>(
            new MathLib::PiecewiseLinearInterpolation(
                std::move(variables), std::move(values), true));
        return std::unique_ptr<RelativePermeability>(
            new ReletivePermeabilityCurve(curve));
    }
    else
    {
        OGS_FATAL(
            "The relative permeability model %s is unavailable.\n"
            "The available models are:"
            "\n\tWettingPhaseVanGenuchten,"
            "\n\tNonWettingPhaseVanGenuchten,"
            "\n\tWettingPhaseBrookCoreyOilGas,"
            "\n\tNonWettingPhaseBrookCoreyOilGas,",
            "\n\tCurve.\n",
            type.data());
    }
}

}  // end namespace
}  // end namespace
