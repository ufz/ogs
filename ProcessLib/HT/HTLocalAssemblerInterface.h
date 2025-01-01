/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on October 11, 2017, 1:35 PM
 */

#pragma once

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{

namespace HT
{
template <typename GlobalDimNodalMatrixType>
struct IntegrationPointData final
{
    IntegrationPointData(GlobalDimNodalMatrixType dNdx_,
                         double const& integration_weight_)
        : dNdx(std::move(dNdx_)), integration_weight(integration_weight_)
    {
    }

    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class HTLocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                  public NumLib::ExtrapolatableElement
{
public:
    HTLocalAssemblerInterface() = default;

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;

    Eigen::Vector3d getFlux(MathLib::Point3d const& pnt_local_coords,
                            double const t,
                            std::vector<double> const& local_x) const override =
        0;
};

}  // namespace HT
}  // namespace ProcessLib
