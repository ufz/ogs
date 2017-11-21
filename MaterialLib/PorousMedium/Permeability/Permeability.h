/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   Permeability.h
 *
 */

#pragma once

#include <Eigen/Dense>

#include "BaseLib/Error.h"
#include "ProcessLib/Parameter/Parameter.h"

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
        ProcessLib::Parameter<double> const& permeability_parameter,
        int const dimension)
        : _permeability_parameter(permeability_parameter), _dimension(dimension)
    {
        if (permeability_parameter.getNumberOfComponents() !=
            _dimension * _dimension)
        {
            OGS_FATAL(
                "The given parameter has %d components, but the permeability "
                "tensor is defined for a %d dimensional problem.",
                permeability_parameter.getNumberOfComponents(), _dimension);
        }
        _intrinsic_permeability_tensor =
            Eigen::MatrixXd(_dimension, _dimension);
    }

    ~Permeability() = default;

    /**
     *  Get property value.
     *  @param t point in time
     *  @param pos spatial position
     *  @param variable    A variable with any double type value.
     *  @param temperature Temperature with any double type value.
     */
    Eigen::MatrixXd const& getValue(
        const double t,
        ProcessLib::SpatialPosition const& pos,
        const double variable,
        const double temperature) const
    {
        (void)variable;
        (void)temperature;

        _intrinsic_permeability_tensor =
            Eigen::Map<Eigen::Matrix<double,
                                     Eigen::Dynamic,
                                     Eigen::Dynamic,
                                     Eigen::RowMajor> const>(
                _permeability_parameter(t, pos).data(), _dimension, _dimension);

        return _intrinsic_permeability_tensor;
    }

private:
    ProcessLib::Parameter<double> const& _permeability_parameter;
    int const _dimension;
    mutable Eigen::MatrixXd _intrinsic_permeability_tensor;
};

}  // end of namespace
}  // end of namespace
