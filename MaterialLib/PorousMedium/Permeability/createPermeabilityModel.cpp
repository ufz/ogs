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

#include "BaseLib/ConfigTree.h"

#include "BaseLib/Error.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
Eigen::MatrixXd createPermeabilityModel(BaseLib::ConfigTree const& config)
{
    auto const values =
        //! \ogs_file_param{material__porous_medium__permeability__values}
        config.getConfigParameter<std::vector<double>>("values");

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

    return MathLib::toMatrix(values, dim, dim);
}

}  // end of namespace
}  // end of namespace
