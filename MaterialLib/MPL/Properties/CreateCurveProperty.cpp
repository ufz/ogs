/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateCurveProperty.h"

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "CurveProperty.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MaterialPropertyLib
{
std::unique_ptr<CurveProperty> createCurveProperty(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    DBUG("Create CurveProperty.");
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Curve");

    //! \ogs_file_param{properties__property__Curve__curve}
    auto curve_name = config.getConfigParameter<std::string>("curve");
    DBUG("Using curve '{:s}'", curve_name.c_str());

    auto const& curve =
        *BaseLib::getOrError(curves, curve_name, "Could not find curve.");

    auto const independent_variable_string =
        //! \ogs_file_param{properties__property__Curve__independent_variable}
        config.getConfigParameter<std::string>("independent_variable");
    DBUG("Using independent_variable '{:s}'",
         independent_variable_string.c_str());
    auto const independent_variable =
        MaterialPropertyLib::convertStringToVariable(
            independent_variable_string);

    return std::make_unique<CurveProperty>(independent_variable, curve);
}

}  // namespace ParameterLib
