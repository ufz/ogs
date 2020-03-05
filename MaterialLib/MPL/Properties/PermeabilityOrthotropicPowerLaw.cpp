/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/PermeabilityOrthotropicPowerLaw.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/CoordinateSystem.h"

namespace MaterialPropertyLib
{
template <int DisplacementDim>
PermeabilityOrthotropicPowerLaw<DisplacementDim>::
    PermeabilityOrthotropicPowerLaw(
        std::array<double, DisplacementDim> intrinsic_permeabilities,
        std::array<double, DisplacementDim>
            exponents,
        ParameterLib::CoordinateSystem const* const local_coordinate_system)
    : _k(std::move(intrinsic_permeabilities)),
      _lambda(std::move(exponents)),
      _local_coordinate_system(local_coordinate_system)
{
}

template <int DisplacementDim>
void PermeabilityOrthotropicPowerLaw<DisplacementDim>::setScale(
    std::variant<Medium*, Phase*, Component*> scale_pointer)
{
    if (std::holds_alternative<Phase*>(scale_pointer))
    {
        _phase = std::get<Phase*>(scale_pointer);
        if (_phase->name != "Solid")
        {
            OGS_FATAL(
                "The property 'PermeabilityOrthotropicPowerLaw' must be "
                "given in the 'Solid' phase, not in '%s' phase.",
                _phase->name.c_str());
        }
    }
    else
    {
        OGS_FATAL(
            "The property 'PermeabilityOrthotropicPowerLaw' is "
            "implemented on the 'phase' scales only.");
    }
}
template <int DisplacementDim>
PropertyDataType PermeabilityOrthotropicPowerLaw<DisplacementDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const /*t*/,
    double const /*dt*/) const
{
    auto const phi =
        std::get<double>(variable_array[static_cast<int>(Variable::porosity)]);
    // TODO (naumov) The phi0 must be evaluated once upon
    // creation/initialization and be stored in a local state.
    // For now assume porosity's initial value does not change with time.
    auto const phi_0 = _phase->property(PropertyType::porosity)
                           .template initialValue<double>(
                               pos, std::numeric_limits<double>::quiet_NaN());

    Eigen::Matrix<double, DisplacementDim, DisplacementDim> k =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();

    Eigen::Matrix<double, DisplacementDim, DisplacementDim> const e =
        _local_coordinate_system == nullptr
            ? Eigen::Matrix<double, DisplacementDim,
                            DisplacementDim>::Identity()
            : _local_coordinate_system->transformation<DisplacementDim>(pos);

    // k = \sum_i k_i (\phi / \phi_0)^{\lambda_i} e_i \otimes e_i
    // e_i \otimes e_i = square matrix e_i,0^2 e_i,0*e_i,1 etc.
    for (int i = 0; i < DisplacementDim; ++i)
    {
        Eigen::Matrix<double, DisplacementDim, DisplacementDim> const
            ei_otimes_ei = e.col(i) * e.col(i).transpose();

        k += _k[i] * std::pow(phi / phi_0, _lambda[i]) * ei_otimes_ei;
    }
    return k;
}
template <int DisplacementDim>
PropertyDataType PermeabilityOrthotropicPowerLaw<DisplacementDim>::dValue(
    VariableArray const& /*variable_array*/, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::strain) &&
           "PermeabilityOrthotropicPowerLaw::dValue is implemented for "
           " derivatives with respect to strain only.");

    return 0;
}

template class PermeabilityOrthotropicPowerLaw<2>;
template class PermeabilityOrthotropicPowerLaw<3>;
}  // namespace MaterialPropertyLib
