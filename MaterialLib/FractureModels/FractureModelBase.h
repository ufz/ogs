/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_FRACTURE_FRACTUREMODELBASE_H_
#define MATERIALLIB_FRACTURE_FRACTUREMODELBASE_H_

#include <Eigen/Eigen>

#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Fracture
{

/**
 * Interface for mechanical fracture material models. Provides updates of the
 * stress for a given current state and also a tangent at that position.
 */
template <int DisplacementDim>
class FractureModelBase
{
public:
    virtual ~FractureModelBase() {}

    /**
     * Computation of the constitutive relation for specific material model.
     * This should be implemented in the derived model.
     *
     * @param t           current time
     * @param x           current position in space
     * @param w_prev      fracture displacement at previous time step
     * @param w           fracture displacement at current time step
     * @param sigma_prev  stress at previous time step
     * @param sigma       stress at current time step
     * @param C           tangent matrix for stress and fracture displacements
     */
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

#endif
