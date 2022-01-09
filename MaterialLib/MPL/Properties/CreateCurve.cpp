/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateCurve.h"

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "Curve.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Curve> createCurve(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Curve");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create Curve {:s}.", property_name);

    //! \ogs_file_param{properties__property__Curve__curve}
    auto curve_name = config.getConfigParameter<std::string>("curve");
    DBUG("Using curve '{:s}'", curve_name);

    auto const& curve =
        *BaseLib::getOrError(curves, curve_name, "Could not find curve.");

    auto const independent_variable_string =
        //! \ogs_file_param{properties__property__Curve__independent_variable}
        config.getConfigParameter<std::string>("independent_variable");
    DBUG("Using independent_variable '{:s}'", independent_variable_string);
    auto const independent_variable =
        MaterialPropertyLib::convertStringToVariable(
            independent_variable_string);

    return std::make_unique<Curve>(
        std::move(property_name), independent_variable, curve);
}

}  // namespace MaterialPropertyLib
