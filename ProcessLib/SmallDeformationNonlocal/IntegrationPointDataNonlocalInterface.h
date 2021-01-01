/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
struct IntegrationPointDataNonlocalInterface;

struct NonlocalIP final
{
    IntegrationPointDataNonlocalInterface* const ip_l_pointer;
    double alpha_kl_times_w_l;
    double distance2;  ///< Squared distance to current integration point.
};

struct IntegrationPointDataNonlocalInterface
{
    virtual ~IntegrationPointDataNonlocalInterface() = default;

    std::vector<NonlocalIP> non_local_assemblers;

    double kappa_d = 0;      ///< damage driving variable.
    double integration_weight;
    double nonlocal_internal_length;
    Eigen::Vector3d coordinates;
    bool active_self = false;
    bool activated = false;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
