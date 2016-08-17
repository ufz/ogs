/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   createPermeabilityModel.cpp
 *
 * Created on August 17, 2016, 2:36 PM
 */

#include "createPermeabilityModel.h"

#include <cassert>

#include "BaseLib/Error.h"
#include "BaseLib/StringTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
CoefMatrix createPermeabilityModel(BaseLib::ConfigTree const* const config)
{
    //! \ogs_file_param{material__porous_medium__permeability__values}
    auto const val_strs = config->getConfigParameter<std::string>("values");
    std::vector<double> values =
        BaseLib::convertStringValues2Values<double>(val_strs);

    auto const data_size = values.size();

    int dim = -1;
    switch (data_size)
    {
        case 1:
            dim = 1;
            break;
        case 4:
            dim = 2;
            break;
        case 9:
            dim = 3;
            break;
        default:
        {
            OGS_FATAL(
                "Number of values for permeability tensor must be 1, 4 or 9.")
        }
    }

    CoefMatrix c_matrix(dim, dim);
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            c_matrix(i, j) = values[i * dim + j];
        }
    }
    return c_matrix;
}

}  // end of namespace
}  // end of namespace
