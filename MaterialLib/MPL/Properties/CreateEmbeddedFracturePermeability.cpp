/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "EmbeddedFracturePermeability.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createEmbeddedFracturePermeability(
    int const geometry_dimension, BaseLib::ConfigTree const& config)
{
    if ((geometry_dimension != 2) && (geometry_dimension != 3))
    {
        OGS_FATAL(
            "The EmbeddedFracturePermeability is implemented only for 2D or 3D "
            "problems");
    }

    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "EmbeddedFracturePermeability");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create EmbeddedFracturePermeability medium property");

    auto const k =
        //! \ogs_file_param{properties__property__EmbeddedFracturePermeability__intrinsic_permeability}
        config.getConfigParameter<double>("intrinsic_permeability");

    auto const b0 =
        //! \ogs_file_param{properties__property__EmbeddedFracturePermeability__initial_aperture}
        config.getConfigParameter<double>("initial_aperture");

    auto const a =
        //! \ogs_file_param{properties__property__EmbeddedFracturePermeability__mean_frac_distance}
        config.getConfigParameter<double>("mean_frac_distance");

    auto const e0 =
        //! \ogs_file_param{properties__property__EmbeddedFracturePermeability__threshold_strain}
        config.getConfigParameter<double>("threshold_strain");

    bool n_const = false;
    Eigen::Matrix<double, 3, 1> n;
    if (auto const n_ptr =
            //! \ogs_file_param{properties__property__EmbeddedFracturePermeability__fracture_normal}
        config.getConfigParameterOptional<std::vector<double>>(
            "fracture_normal"))
    {
        if ((*n_ptr).size() != 3)
        {
            OGS_FATAL(
                "The size of the fracture normal vector must be 3, but is %d.",
                (*n_ptr).size());
        }
        DBUG("Using constant fracture normal vector.");
        std::copy_n((*n_ptr).data(), 3, n.data());
        n_const = true;
        n /= n.norm();
    }
    else
    {
        DBUG(
            "No constant fracture normal was given. By default it will be "
            "determined as the third principal stress vector.");
    }

    auto const jf =
        //! \ogs_file_param{properties__property__EmbeddedFracturePermeability__jacobian_factor}
        config.getConfigParameter<double>("jacobian_factor", 0.);

    if (geometry_dimension == 2)
    {
        return std::make_unique<EmbeddedFracturePermeability<2>>(
            std::move(property_name), n, n_const, k, b0, a, e0, jf);
    }
    return std::make_unique<EmbeddedFracturePermeability<3>>(
        std::move(property_name), n, n_const, k, b0, a, e0, jf);
}
}  // namespace MaterialPropertyLib
