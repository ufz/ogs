// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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
