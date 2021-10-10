/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSimplifiedElasticityModel.h"

#include <Eigen/Dense>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "HydrostaticElasticityModel.h"
#include "RigidElasticityModel.h"
#include "SimplifiedElasticityModel.h"
#include "UniaxialElasticityModel.h"
#include "UserDefinedElasticityModel.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
std::unique_ptr<SimplifiedElasticityModel> createElasticityModel(
    BaseLib::ConfigTree const& config)
{
    std::unique_ptr<SimplifiedElasticityModel> simplified_elasticity =
        std::make_unique<RigidElasticityModel>();
    if (auto const simplified_elasticity_switch =
            //! \ogs_file_param{prj__processes__process__THERMO_RICHARDS_FLOW__simplified_elasticity}
        config.getConfigParameterOptional<std::string>("simplified_elasticity"))
    {
        DBUG("Using simplified_elasticity for the Richards flow equation");
        if (*simplified_elasticity_switch == "uniaxial")
        {
            DBUG("assuming local uniaxial deformation only.");
            simplified_elasticity = std::make_unique<UniaxialElasticityModel>();
        }
        else if (*simplified_elasticity_switch == "hydrostatic")
        {
            DBUG("assuming constant hydrostatic stress locally.");
            simplified_elasticity =
                std::make_unique<HydrostaticElasticityModel>();
        }
        else if (*simplified_elasticity_switch == "user_defined")
        {
            DBUG("using user defined elasticity model.");
            simplified_elasticity =
                std::make_unique<UserDefinedElasticityModel>();
        }
        else if (*simplified_elasticity_switch == "rigid")
        {
            DBUG("using user defined elasticity model.");
            simplified_elasticity = std::make_unique<RigidElasticityModel>();
        }
    }
    return simplified_elasticity;
}
}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
