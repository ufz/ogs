/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_EXTRAPOLATION_EXTRAPOLATABLEELEMENT_H
#define NUMLIB_EXTRAPOLATION_EXTRAPOLATABLEELEMENT_H

#include <Eigen/Core>

namespace NumLib
{
/*! Interface for providing shape matrices and integration point values for
 *  extrapolation,
 *
 * Local assemblers that want to have some integration point values extrapolated
 * using an Extrapolator have to implement this interface.
 */
class ExtrapolatableElement
{
public:
    //! Provides the shape matrix at the given integration point.
    virtual Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const = 0;

    virtual ~ExtrapolatableElement() = default;
};

}  // namespace NumLib

#endif  // NUMLIB_EXTRAPOLATION_EXTRAPOLATABLEELEMENT_H
