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
enum class BrookCoreyOrVanGenuchtenModelType
{
    WettingPhaseVanGenuchten,
    NonWettingPhaseVanGenuchten,
    WettingPhaseBrookCoreyOilGas,
    NonWettingPhaseBrookCoreyOilGas
};

/**
 *     \param config ConfigTree object which contains the input data
 *                   including <type>WettingPhaseBrookCoreyOilGas</type>
 *                          or <type>NonWettingPhaseBrookCoreyOilGas</type>
 *                          or <type>WettingPhaseVanGenuchten</type>
 *                          or <type>NonWettingPhaseVanGenuchten</type>
 *                   and it has a tag of <relative_permeability>
*/
static std::unique_ptr<RelativePermeability> createBrookCoreyOrVanGenuchten(
    BaseLib::ConfigTree const& config,
    const BrookCoreyOrVanGenuchtenModelType mode_type)
{
    std::array<double, 4> parameters = {
        {//! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type__sr}
         config.getConfigParameter<double>("sr"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type__smax}
         config.getConfigParameter<double>("smax"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type__m}
         config.getConfigParameter<double>("m"),
         //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type__krel_min}
         config.getConfigParameter<double>("krel_min")}};

    switch (mode_type)
    {
        case BrookCoreyOrVanGenuchtenModelType::WettingPhaseVanGenuchten:
            assert(parameters[2] <= 1.0);  // m <= 1
            return std::unique_ptr<RelativePermeability>(
                new WettingPhaseVanGenuchten(parameters));
        case BrookCoreyOrVanGenuchtenModelType::NonWettingPhaseVanGenuchten:
            assert(parameters[2] <= 1.0);  // m <= 1
            return std::unique_ptr<RelativePermeability>(
                new NonWettingPhaseVanGenuchten(parameters));
        case BrookCoreyOrVanGenuchtenModelType::WettingPhaseBrookCoreyOilGas:
            assert(parameters[2] >= 1.0);  // m >= 1
            return std::unique_ptr<RelativePermeability>(
                new WettingPhaseBrookCoreyOilGas(parameters));
        case BrookCoreyOrVanGenuchtenModelType::NonWettingPhaseBrookCoreyOilGas:
            assert(parameters[2] >= 1.0);  // m >= 1
            return std::unique_ptr<RelativePermeability>(
                new NonWettingPhaseBrookCoreyOilGas(parameters));
    }

    return nullptr;
}

std::unique_ptr<RelativePermeability> createRelativePermeabilityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "WettingPhaseVanGenuchten")
    {
        return createBrookCoreyOrVanGenuchten(
            config,
            BrookCoreyOrVanGenuchtenModelType::WettingPhaseVanGenuchten);
    }
    else if (type == "NonWettingPhaseVanGenuchten")
    {
        return createBrookCoreyOrVanGenuchten(
            config,
            BrookCoreyOrVanGenuchtenModelType::NonWettingPhaseVanGenuchten);
    }
    else if (type == "WettingPhaseBrookCoreyOilGas")
    {
        return createBrookCoreyOrVanGenuchten(
            config,
            BrookCoreyOrVanGenuchtenModelType::WettingPhaseBrookCoreyOilGas);
    }
    else if (type == "NonWettingPhaseBrookCoreyOilGas")
    {
        return createBrookCoreyOrVanGenuchten(
            config,
            BrookCoreyOrVanGenuchtenModelType::NonWettingPhaseBrookCoreyOilGas);
    }
    else if (type == "Curve")
    {
        std::vector<double> variables, values;
        //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type__curve}
        auto const& curve_config = config.getConfigSubtree("curve");
        //! \ogs_file_param{material_property__porous_medium__porous_medium__relative_permeability__type__curve__data}
        for (auto const& data_string :
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
            "The capillary pressure model %s is unavailable.\n"
            "The available models are:"
            "\n\tWettingPhaseVanGenuchten,"
            "\n\tNonWettingPhaseVanGenuchten,"
            "\n\tWettingPhaseBrookCoreyOilGas,"
            "\n\tNonWettingPhaseBrookCoreyOilGas.\n",
            type.data());
    }
}

}  // end namespace
}  // end namespace
