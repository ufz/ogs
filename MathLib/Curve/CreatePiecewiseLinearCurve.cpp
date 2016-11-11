/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreatePiecewiseLinearCurve.cpp
 *
 * Created on November 11, 2016, 10:49 AM
 */

#include "CreatePiecewiseLinearCurve.h"

#include "BaseLib/Error.h"
#include "BaseLib/ConfigTree.h"

#include "PiecewiseLinearCurve.h"

namespace MathLib
{
std::unique_ptr<PiecewiseLinearCurve> createPiecewiseLinearCurve(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{curve__type}
    config.checkConfigParameter("type", "PiecewiseLinear");
    std::vector<double> x, y;
    for (auto const& point_config :
         //! \ogs_file_param{curve__point}
         config.getConfigSubtreeList("point"))
    {
        const auto& point =
            //! \ogs_file_param{curve__point__data}
            point_config.getConfigParameter<std::vector<double>>("data");
        assert(point.size() == 2);
        x.push_back(point[0]);
        y.push_back(point[1]);
    }

    return std::unique_ptr<MathLib::PiecewiseLinearCurve>(
        new MathLib::PiecewiseLinearCurve(std::move(x), std::move(y)));
}
}
