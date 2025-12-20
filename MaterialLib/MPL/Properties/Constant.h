// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/// The constant property class. This property simply retrieves the stored
/// constant value. It accepts all datatypes defined in PropertyDataType
/// (currently: double, Vector, Tensor, std::string)
class Constant final : public Property
{
public:
    /// This constructor accepts single values of any data type defined in the
    /// PropertyDataType definition and sets the protected attribute value_ of
    /// the base class Property to that value.
    Constant(std::string name, PropertyDataType const& v);

    /// If the property is no a scalar constant, the returned scalar zero value
    /// is used to form a tensor or a vector via calling formEigenTensor or
    /// FormEigenVector, respectively
    PropertyDataType dValue(VariableArray const& /*variable_array*/,
                            VariableArray const& /*variable_array_prev*/,
                            Variable const /*variable*/,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override
    {
        return 0.0;
    }
};
}  // namespace MaterialPropertyLib
