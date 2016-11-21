/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CurveScaledParameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createCurveScaledParameter(
    std::string const& name,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "CurveScaled");

    //! \ogs_file_param{prj__parameters__parameter__CurveScaled__curve}
    auto curve_name = config.getConfigParameter<std::string>("curve");
    DBUG("Using curve %s", curve_name.c_str());

    auto const curve_it = curves.find(curve_name);
    if (curve_it == curves.end())
        OGS_FATAL("Curve `%s' does not exists.", curve_name.c_str());

    //! \ogs_file_param{prj__parameters__parameter__CurveScaled__parameter}
    auto referenced_parameter_name =
            config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", referenced_parameter_name.c_str());

    // TODO other data types than only double
    return std::unique_ptr<ParameterBase>(new CurveScaledParameter<double>(
        name, *curve_it->second, referenced_parameter_name));
}

}  // ProcessLib
