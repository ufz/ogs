/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Dense>

#include "BaseLib/Error.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace PorousMedium
{
/// The class implements a basic permeability model that employs a parameter
/// (for instance a constant parameter or mesh cell dependend parameter) to fill
/// the intrinsic permeability tensor.
class Permeability
{
public:
    explicit Permeability(
        ParameterLib::Parameter<double> const& permeability_parameter,
        int const dimension)
        : permeability_parameter_(permeability_parameter), dimension_(dimension)
    {
        if (permeability_parameter.getNumberOfComponents() !=
            dimension_ * dimension_)
        {
            OGS_FATAL(
                "The given parameter has {:d} components, but the permeability "
                "tensor is defined for a {:d} dimensional problem.",
                permeability_parameter.getNumberOfComponents(), dimension_);
        }
    }

    virtual ~Permeability() = default;

    /**
     *  Get the intrinsic permeability tensor.
     *  @param t point in time
     *  @param pos spatial position
     *  @param variable    A variable with any double type value.
     *  @param temperature Temperature with any double type value.
     */
    virtual Eigen::MatrixXd getValue(const double t,
                                     ParameterLib::SpatialPosition const& pos,
                                     const double variable,
                                     const double temperature) const
    {
        (void)variable;
        (void)temperature;

        return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                        Eigen::RowMajor> const>(
            permeability_parameter_(t, pos).data(), dimension_, dimension_);
    }

private:
    ParameterLib::Parameter<double> const& permeability_parameter_;
    int const dimension_;
};

}  // namespace PorousMedium
}  // namespace MaterialLib
