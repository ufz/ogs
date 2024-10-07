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

    std::ifstream file(path_to_file, std::ios::binary);
    if (!file)
    {
        OGS_FATAL("File {:s} at path {:s} for curve definition not found",
                  filename, prj_dir);
    }

    // Determine file size
    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    // Initialize vector with the right size
    std::vector<double> result(size / sizeof(double));

    // read data directly into the vector
    if (!file.read(reinterpret_cast<char*>(result.data()), size))
    {
        OGS_FATAL("Could not read curve definition from file {:s}.", filename);
    }

    if constexpr (std::endian::native != std::endian::little)
    {
        // swap endians if needed
        std::transform(
            result.begin(), result.end(), result.begin(),
            [](double value)
            {
                return std::bit_cast<double>(boost::endian::endian_reverse(
                    std::bit_cast<uint64_t>(value)));
            });
    }

    return result;
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
