/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   createPermeabilityModel.h
 *
 * Created on August 17, 2016, 2:36 PM
 */

#ifndef CREATEPERMEABILITYMODEL_H
#define CREATEPERMEABILITYMODEL_H

#include <Eigen/Dense>

#include "BaseLib/ConfigTree.h"

namespace MaterialLib
{
namespace PorousMedium
{
using CoefMatrix = Eigen::MatrixXd;

/// Create a porosity model
/// \param config  ConfigTree object has a tag of <permeability>
CoefMatrix createPermeabilityModel(BaseLib::ConfigTree const& config);

}  // end of namespace
}  // end of namespace

#endif /* CREATEPERMEABILITYMODEL_H */
