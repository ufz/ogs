/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
        : _permeability_parameter(permeability_parameter), _dimension(dimension)
    {
        if (permeability_parameter.getNumberOfGlobalComponents() !=
            _dimension * _dimension)
        {
            OGS_FATAL(
                "The given parameter has {:d} components, but the permeability "
                "tensor is defined for a {:d} dimensional problem.",
                permeability_parameter.getNumberOfGlobalComponents(),
                _dimension);
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
            _permeability_parameter(t, pos).data(), _dimension, _dimension);
    }

private:
    ParameterLib::Parameter<double> const& _permeability_parameter;
    int const _dimension;
};

}  // namespace PorousMedium
}  // namespace MaterialLib
