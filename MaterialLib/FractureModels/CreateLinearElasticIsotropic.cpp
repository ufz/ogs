/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateLinearElasticIsotropic.h"

#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter
#include "LinearElasticIsotropic.h"

namespace MaterialLib
{
namespace Fracture
{

template <int DisplacementDim>
std::unique_ptr<FractureModelBase<DisplacementDim>>
createLinearElasticIsotropic(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "LinearElasticIsotropic");
    DBUG("Create LinearElasticIsotropic material");

    auto& Kn = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticIsotropic__normal_stiffness}
        config, "normal_stiffness", parameters, 1);

    auto& Ks = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticIsotropic__shear_stiffness}
        config, "shear_stiffness", parameters, 1);

    typename LinearElasticIsotropic<DisplacementDim>::MaterialProperties mp{
        Kn, Ks};

    return std::unique_ptr<LinearElasticIsotropic<DisplacementDim>>{
        new LinearElasticIsotropic<DisplacementDim>{mp}};
}


template
std::unique_ptr<FractureModelBase<2>>
createLinearElasticIsotropic(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

template
std::unique_ptr<FractureModelBase<3>>
createLinearElasticIsotropic(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}  // namespace Fracture
}  // namespace MaterialLib
