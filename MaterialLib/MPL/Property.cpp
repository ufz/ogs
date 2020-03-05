/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Property.h"

#include <string>

namespace MaterialPropertyLib
{
PropertyDataType fromVector(std::vector<double> const& values)
{
    switch (values.size())
    {
        case 1:
        {
            return values[0];
        }
        case 2:
        {
            return Eigen::Vector2d{values[0], values[1]};
        }
        case 3:
        {
            return Eigen::Vector3d{values[0], values[1], values[2]};
        }
        case 4:
        {
            using M = Eigen::Matrix2d;
            return M{Eigen::Map<M const>{values.data(), 2, 2}};
        }
        case 6:
        {
            // Symmetric Tensor - xx, yy, zz, xy, xz, yz
            using M = Eigen::Matrix<double, 6, 1>;
            return M{Eigen::Map<M const>{values.data(), 6}};
        }
        case 9:
        {
            using M = Eigen::Matrix3d;
            return M{Eigen::Map<M const>{values.data(), 3, 3}};
        }
        default:
        {
            OGS_FATAL(
                "Conversion of a %d-vector to PropertyDataType is not "
                "implemented.",
                values.size());
        }
    }
}

PropertyDataType Property::initialValue(
    ParameterLib::SpatialPosition const& pos, double const t) const
{
    return value(VariableArray{}, pos, t,
                 std::numeric_limits<double>::quiet_NaN());
}

PropertyDataType Property::value() const
{
    return _value;
}

/// The default implementation of this method only returns the property value
/// without altering it.
PropertyDataType Property::value(VariableArray const& /*variable_array*/,
                                 ParameterLib::SpatialPosition const& /*pos*/,
                                 double const /*t*/, double const /*dt*/) const
{
    return _value;
}

/// The default implementation of this method only returns the
/// property value derivative without altering it.
PropertyDataType Property::dValue(VariableArray const& /*variable_array*/,
                                  Variable const /*variable*/,
                                  ParameterLib::SpatialPosition const& /*pos*/,
                                  double const /*t*/, double const /*dt*/) const
{
    return _dvalue;
}

/// Default implementation: 2nd derivative of any constant property is zero.
PropertyDataType Property::d2Value(VariableArray const& /*variable_array*/,
                                   Variable const /*variable*/,
                                   Variable const /*variable*/,
                                   ParameterLib::SpatialPosition const& /*pos*/,
                                   double const /*t*/,
                                   double const /*dt*/) const
{
    return 0.0;
}
}  // namespace MaterialPropertyLib
