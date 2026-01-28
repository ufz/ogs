// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreatePiecewiseLinearCurve.h"

#include <boost/endian/conversion.hpp>
#include <fstream>
#include <iostream>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

namespace MathLib
{
PiecewiseLinearCurveConfig parsePiecewiseLinearCurveConfig(
    BaseLib::ConfigTree const& config)
{
    const bool read_from_file =  //! \ogs_file_param{curve__read_from_file}
        config.getConfigParameter<bool>("read_from_file", false);

    std::vector<double> x;
    std::vector<double> y;

    if (read_from_file == true)
    {
        auto const coords_file_name =
            //! \ogs_file_param{curve__coords}
            config.getConfigParameter<std::string>("coords");

        auto const values_file_name =
            //! \ogs_file_param{curve__values}
            config.getConfigParameter<std::string>("values");

        x = BaseLib::readDoublesFromBinaryFile(
            coords_file_name, config.projectDirectory().string());

        y = BaseLib::readDoublesFromBinaryFile(
            values_file_name, config.projectDirectory().string());
    }
    else
    {
        x =
            //! \ogs_file_param{curve__coords}
            config.getConfigParameter<std::vector<double>>("coords");
        y =
            //! \ogs_file_param{curve__values}
            config.getConfigParameter<std::vector<double>>("values");
    }

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

    return {std::move(x), std::move(y)};
}
}  // namespace MathLib
