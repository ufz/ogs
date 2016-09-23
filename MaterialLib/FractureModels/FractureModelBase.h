/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Fracture
{

template <int DisplacementDim>
class FractureModelBase
{
public:
    virtual ~FractureModelBase() {}

    virtual void computeConstitutiveRelation(
            double const t,
            ProcessLib::SpatialPosition const& x,
            Eigen::Ref<Eigen::VectorXd const> w_prev,
            Eigen::Ref<Eigen::VectorXd const> w,
            Eigen::Ref<Eigen::VectorXd const> sigma_prev,
            Eigen::Ref<Eigen::VectorXd> sigma,
            Eigen::Ref<Eigen::MatrixXd> C)  = 0;

};

}  // namespace Fracture
}  // namespace MaterialLib

