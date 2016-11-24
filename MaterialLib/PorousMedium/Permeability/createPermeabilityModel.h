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

#ifndef OGS_CREATEPERMEABILITYMODEL_H
#define OGS_CREATEPERMEABILITYMODEL_H

#include <Eigen/Dense>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace PorousMedium
{
/** Create a porosity model
 *  @param config  ConfigTree object has a tag of `<permeability>`
 */
Eigen::MatrixXd createPermeabilityModel(BaseLib::ConfigTree const& config);

}  // end of namespace
}  // end of namespace

#endif /* CREATEPERMEABILITYMODEL_H */
