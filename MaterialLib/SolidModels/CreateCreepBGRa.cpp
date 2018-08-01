/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateCreepBGRa.cpp
 *  Created on July 11, 2018, 2:26 PM
 */

#include "CreateCreepBGRa.h"

#include "CreateLinearElasticIsotropic.h"
#include "CreateNewtonRaphsonSolverParameters.h"

#include "CreepBGRa.h"

#include "MechanicsBase.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
namespace Creep
{
template <int DisplacementDim>
std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
createCreepBGRa(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "CreepBGRa");
    DBUG("Create CreepBGRa material");

    // Read elastic data frist.
    const bool skip_type_checking = true;
    auto elastic_data =
        MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
            parameters, config, skip_type_checking);

    //! \ogs_file_param_special{material__solid__constitutive_relation__CreepBGRa__a}
    const auto A = config.getConfigParameter<double>("a");
    DBUG("CreepBGRa parameter A=%g", A);

    //! \ogs_file_param_special{material__solid__constitutive_relation__CreepBGRa__n}
    const auto n = config.getConfigParameter<double>("n");
    DBUG("CreepBGRa parameter n=%g", n);

    //! \ogs_file_param_special{material__solid__constitutive_relation__CreepBGRa__sigma0}
    const auto sigma0 = config.getConfigParameter<double>("sigma0");
    DBUG("CreepBGRa parameter sigma0=%g", sigma0);

    //! \ogs_file_param_special{material__solid__constitutive_relation__CreepBGRa__q}
    const auto Q = config.getConfigParameter<double>("q");
    DBUG("CreepBGRa parameter Q=%g", Q);

    auto const nonlinear_solver_parameters =
        createNewtonRaphsonSolverParameters(config);

    return std::unique_ptr<CreepBGRa<DisplacementDim>>{
        new CreepBGRa<DisplacementDim>{elastic_data->getMaterialProperties(),
                                       nonlinear_solver_parameters, A, n,
                                       sigma0, Q}};
}

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>
createCreepBGRa<2>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>
createCreepBGRa<3>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);
}  // namespace Creep
}  // namespace Solids
}  // namespace MaterialLib
