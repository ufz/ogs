/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 6, 2023, 3:38 PM
 */

#include "CreateLinearElasticTransverseIsotropic.h"

#include "ParameterLib/Utils.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<LinearElasticTransverseIsotropic<DisplacementDim>>
createLinearElasticTransverseIsotropic(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config, const bool skip_type_checking)
{
    if (!skip_type_checking)
    {
        //! \ogs_file_param{material__solid__constitutive_relation__type}
        config.checkConfigParameter("type", "LinearElasticTransverseIsotropic");
        DBUG("Create LinearElasticTransverseIsotropic material");
    }

    auto const& E_i = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticTransverseIsotropic__youngs_modulus_i}
        config, "youngs_modulus_i", parameters, 1);
    DBUG("Use '{}' as the in-plane Young’s modulus, E_i.", E_i.name);

    auto const& E_a = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticTransverseIsotropic__youngs_modulus_a}
        config, "youngs_modulus_a", parameters, 1);
    DBUG(
        "Use '{}' as the Young’s modulus w.r.t. the direction of anisotropy, "
        "E_a.",
        E_a.name);

    auto const& nu_ii = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticTransverseIsotropic__poissons_ratio_ii}
        config, "poissons_ratio_ii", parameters, 1);
    DBUG("Use '{}' as the in-plane Poisson’s ratio, nu_ii.", nu_ii.name);

    auto const& nu_ia = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticTransverseIsotropic__poissons_ratio_ia}
        config, "poissons_ratio_ia", parameters, 1);
    DBUG(
        "Use '{}' as the Poisson ratio perpendicular to the plane of isotropy, "
        "due to strain in the plane of isotropy, nu_ia.",
        nu_ia.name);

    auto const& G_ia = ParameterLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__LinearElasticTransverseIsotropic__shear_modulus_ia}
        config, "shear_modulus_ia", parameters, 1);
    DBUG(
        "Use '{}' as the shear modulus between directions of isotropy and "
        "anisotropy , G_ia.",
        G_ia.name);

    return std::make_unique<LinearElasticTransverseIsotropic<DisplacementDim>>(
        E_i, E_a, nu_ii, nu_ia, G_ia, local_coordinate_system);
}

template std::unique_ptr<LinearElasticTransverseIsotropic<2>>
createLinearElasticTransverseIsotropic(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config, const bool skip_type_checking);

template std::unique_ptr<LinearElasticTransverseIsotropic<3>>
createLinearElasticTransverseIsotropic(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config, const bool skip_type_checking);
}  // namespace Solids
}  // namespace MaterialLib
