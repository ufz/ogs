/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "OrthotropicEmbeddedFracturePermeability.h"
#include "ParameterLib/Utils.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createOrthotropicEmbeddedFracturePermeability(
    int const geometry_dimension, BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    if ((geometry_dimension != 2) && (geometry_dimension != 3))
    {
        OGS_FATAL(
            "The OrthotropicEmbeddedFracturePermeability is implemented only "
            "for 2D or 3D problems");
    }

    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "OrthotropicEmbeddedFracturePermeability");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create OrthotropicEmbeddedFracturePermeability medium property");

    auto const a_i =
        //! \ogs_file_param{properties__property__OrthotropicEmbeddedFracturePermeability__mean_frac_distances}
        config.getConfigParameter<std::vector<double>>("mean_frac_distances");
    if (a_i.size() != 3)
    {
        OGS_FATAL(
            "The size of the mean fracture distances vector must be 3, but is "
            "{}.",
            a_i.size());
    }

    auto const e_i0 =
        //! \ogs_file_param{properties__property__OrthotropicEmbeddedFracturePermeability__threshold_strains}
        config.getConfigParameter<std::vector<double>>("threshold_strains");
    if (e_i0.size() != 3)
    {
        OGS_FATAL(
            "The size of the mean threshold strains vector must be 3, but is "
            "{}.",
            e_i0.size());
    }

    auto const n =
        //! \ogs_file_param{properties__property__OrthotropicEmbeddedFracturePermeability__fracture_normals}
        config.getConfigParameter<std::vector<double>>("fracture_normals");
    if (n.size() != 6)
    {
        OGS_FATAL(
            "The size of the fracture normals vector must be 6, but is {}.",
            n.size());
    }
    Eigen::Vector3d const n1 = Eigen::Vector3d({n[0], n[1], n[2]}).normalized();
    Eigen::Vector3d const n2 = Eigen::Vector3d({n[3], n[4], n[5]}).normalized();

    if (n1.dot(n2) > std::numeric_limits<double>::epsilon())
    {
        OGS_FATAL(
            "The given fracture normals are not orthogonal. Please provide two "
            "orthogonal fracture normals");
    }

    Eigen::Matrix3d const n_i =
        (Eigen::Matrix3d() << n1, n2, n1.cross(n2)).finished();

    std::string const intrinsic_permeability_param_name =
        //! \ogs_file_param{properties__property__OrthotropicEmbeddedFracturePermeability__intrinsic_permeability}
        config.getConfigParameter<std::string>("intrinsic_permeability");

    auto const& k = ParameterLib::findParameter<double>(
        intrinsic_permeability_param_name, parameters, 0, nullptr);

    std::string const fracture_rotation_xy_param_name =
        //! \ogs_file_param{properties__property__OrthotropicEmbeddedFracturePermeability__fracture_rotation_xy}
        config.getConfigParameter<std::string>("fracture_rotation_xy");

    auto const& phi_xy = ParameterLib::findParameter<double>(
        fracture_rotation_xy_param_name, parameters, 0, nullptr);

    std::string const fracture_rotation_yz_param_name =
        //! \ogs_file_param{properties__property__OrthotropicEmbeddedFracturePermeability__fracture_rotation_yz}
        config.getConfigParameter<std::string>("fracture_rotation_yz");

    auto const& phi_yz = ParameterLib::findParameter<double>(
        fracture_rotation_yz_param_name, parameters, 0, nullptr);

    auto const jf =
        //! \ogs_file_param{properties__property__OrthotropicEmbeddedFracturePermeability__jacobian_factor}
        config.getConfigParameter<double>("jacobian_factor", 0.);

    if (geometry_dimension == 2)
    {
        return std::make_unique<OrthotropicEmbeddedFracturePermeability<2>>(
            std::move(property_name), a_i, e_i0, n_i, k, phi_xy, phi_yz, jf);
    }
    return std::make_unique<OrthotropicEmbeddedFracturePermeability<3>>(
        std::move(property_name), a_i, e_i0, n_i, k, phi_xy, phi_yz, jf);
}
}  // namespace MaterialPropertyLib
