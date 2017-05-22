/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreatePiecewiseLinearCurve-impl.h
 *
 * Created on November 11, 2016, 10:49 AM
 */

#include "CreatePiecewiseLinearCurve.h"

#include "BaseLib/Error.h"
#include "BaseLib/ConfigTree.h"

namespace MathLib
{
template <typename CurveType>
std::unique_ptr<CurveType> createPiecewiseLinearCurve(
    BaseLib::ConfigTree const& config)
{
    auto x =
        //! \ogs_file_param{curve__coords}
        config.getConfigParameter<std::vector<double>>("coords");
    auto y =
        //! \ogs_file_param{curve__values}
        config.getConfigParameter<std::vector<double>>("values");

    if (x.empty() || y.empty())
    {
        OGS_FATAL("The given coordinates or values vector is empty.");
    }
    if (x.size() != y.size())
    {
        OGS_FATAL(
            "The given coordinates and values vector sizes are "
            "different.");
    }
    return std::make_unique<CurveType>(std::move(x), std::move(y));
}
}
