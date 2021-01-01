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
namespace Deformation
{
/// Divergence of displacement, the volumetric strain.
template <int DisplacementDim, int NPOINTS, typename DNDX_Type>
double divergence(
    const Eigen::Ref<Eigen::Matrix<double, NPOINTS * DisplacementDim, 1> const>&
        u,
    DNDX_Type const& dNdx)
{
    double divergence = 0;
    for (int i = 0; i < DisplacementDim; ++i)
    {
        divergence += dNdx.template block<1, NPOINTS>(i, 0) *
                      u.template segment<NPOINTS>(i * NPOINTS);
    }
    return divergence;
}
}  // namespace Deformation
}  // namespace ProcessLib
