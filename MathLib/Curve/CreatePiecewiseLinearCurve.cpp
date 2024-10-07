/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on November 11, 2016, 10:49 AM
 */

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

std::vector<double> readDoublesFromBinaryFile(const std::string& filename)
{
    auto prj_dir = BaseLib::getProjectDirectory();
    std::string path_to_file = BaseLib::joinPaths(prj_dir, filename);
    std::string file_extension = BaseLib::getFileExtension(filename);
    if (file_extension != ".bin")
    {
        OGS_FATAL(
            "Currently only binary files with extension '.bin' supported. The "
            "specified file has extension {:s}.",
            file_extension)
    }
    return BaseLib::readBinaryVector<double>(path_to_file);
}

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

        x = readDoublesFromBinaryFile(coords_file_name);

        y = readDoublesFromBinaryFile(values_file_name);
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
