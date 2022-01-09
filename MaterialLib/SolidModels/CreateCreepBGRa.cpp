/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on July 11, 2018, 2:26 PM
 */

#include "CreateCreepBGRa.h"

#include "BaseLib/ConfigTree.h"
#include "CreateLinearElasticIsotropic.h"
#include "CreepBGRa.h"
#include "MechanicsBase.h"
#include "NumLib/CreateNewtonRaphsonSolverParameters.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"

namespace MaterialLib
{
namespace Solids
{
namespace Creep
{
template <int DisplacementDim>
std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
createCreepBGRa(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "CreepBGRa");
    DBUG("Create CreepBGRa material");

    // Read elastic data first.
    const bool skip_type_checking = true;
    auto elastic_data =
        MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
            parameters, config, skip_type_checking);

    auto& A = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__CreepBGRa__a}
        config, "a", parameters, 1);

    auto& n = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__CreepBGRa__n}
        config, "n", parameters, 1);

    auto& sigma0 = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__CreepBGRa__sigma0}
        config, "sigma0", parameters, 1);

    auto& Q = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__CreepBGRa__q}
        config, "q", parameters, 1);

    auto const& nonlinear_solver_config =
        //! \ogs_file_param{material__solid__constitutive_relation__CreepBGRa__nonlinear_solver}
        config.getConfigSubtree("nonlinear_solver");
    auto const nonlinear_solver_parameters =
        NumLib::createNewtonRaphsonSolverParameters(nonlinear_solver_config);

    return std::unique_ptr<CreepBGRa<DisplacementDim>>{
        new CreepBGRa<DisplacementDim>{elastic_data->getMaterialProperties(),
                                       nonlinear_solver_parameters, A, n,
                                       sigma0, Q}};
}

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>
createCreepBGRa<2>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>
createCreepBGRa<3>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);
}  // namespace Creep
}  // namespace Solids
}  // namespace MaterialLib
