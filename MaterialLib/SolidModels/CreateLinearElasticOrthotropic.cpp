/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearElasticOrthotropic.h"
#include "ParameterLib/Utils.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<LinearElasticOrthotropic<DisplacementDim>>
createLinearElasticOrthotropic(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config, const bool skip_type_checking)
{
    if (!skip_type_checking)
    {
        //! \ogs_file_param{material__solid__constitutive_relation__type}
        config.checkConfigParameter("type", "LinearElasticOrthotropic");
        DBUG("Create LinearElasticOrthotropic material");
    }

    // The three Youngs moduli are required even in two-dimensional case. Same
    // for the shear moduli and the Poissons ratios.
    auto& youngs_moduli = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticOrthotropic__youngs_moduli}
        config, "youngs_moduli", parameters, 3);
    DBUG("Use '{:s}' as youngs_moduli parameter.", youngs_moduli.name);

    // Shear moduli
    auto& shear_moduli = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticOrthotropic__shear_moduli}
        config, "shear_moduli", parameters, 3);
    DBUG("Use '{:s}' as shear_moduli parameter.", shear_moduli.name);

    // Poissons ratios
    auto& poissons_ratios = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticOrthotropic__poissons_ratios}
        config, "poissons_ratios", parameters, 3);
    DBUG("Use '{:s}' as poissons_ratios parameter.", poissons_ratios.name);

    typename LinearElasticOrthotropic<DisplacementDim>::MaterialProperties mp{
        youngs_moduli, shear_moduli, poissons_ratios};

    return std::make_unique<LinearElasticOrthotropic<DisplacementDim>>(
        mp, local_coordinate_system);
}

template std::unique_ptr<LinearElasticOrthotropic<2>>
createLinearElasticOrthotropic(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config, const bool skip_type_checking);

template std::unique_ptr<LinearElasticOrthotropic<3>>
createLinearElasticOrthotropic(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config, const bool skip_type_checking);
}  // namespace Solids
}  // namespace MaterialLib
