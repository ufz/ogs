/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 12, 2021, x:xx AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

namespace ParameterLib
{
struct CoordinateSystem;
template <typename T>
struct Parameter;
}  // namespace ParameterLib

namespace MaterialPropertyLib
{
/**
 *  \brief A gas pressure dependent intrinsic permeability model.
 *
 *  The model was proposed
 *  in \cite xu2011simulation and it was further investigated
 *  in \cite xu2013coupled .
 *
 *  The model takes the form of
 *  \f[\mathbf{k} = f(p_g)\, \mathbf{k}_0
 *  =
 *  \left\{
 *  \begin{array}{cc}
 *  (1+a_1 p_g) \mathbf{k}_0 & p_g \leq p_\text{thr} \\
 *  (a_2 (p_g - p_\text{thr}) + 1 + a_1 p_\text{thr}) \mathbf{k}_0 & p_g >
 * p_\text{thr} \end{array} \right.\f]
 *  where \f$\mathbf{k}\f$ is the permeability, \f$\mathbf{k}_0\f$ is the
 * initial intrinsic permeability, \f$p_g\f$ is the gas pressure, \f$a_1\f$,
 * \f$a_2\f$ and \f$p_\text{thr}\f$ are three parameters.
 *
 */
template <int DisplacementDim>
class GasPressureDependentPermeability final : public Property
{
public:
    GasPressureDependentPermeability(
        std::string name, ParameterLib::Parameter<double> const& k0,
        double const a1, double const a2, double const pressure_threshold,
        double const minimum_permeability, double const maximum_permeability,
        ParameterLib::CoordinateSystem const* const local_coordinate_system);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    /// Initial intrinsic permeability.
    ParameterLib::Parameter<double> const& k0_;
    /// A factor taking into account the influence of the gas pressure on the
    /// permeability.
    double const a1_;
    /// A factor taking into account the influence of the gas overpressure
    /// (difference of the gas pressure and the threshold pressure) on the
    /// permeability.
    double const a2_;
    /// The threshold pressure which defines the low and high gas pressure
    /// domain. Different permeability functions are defined for each domain.
    double const pressure_threshold_;
    double const minimum_permeability_;
    double const maximum_permeability_;
    ParameterLib::CoordinateSystem const* const local_coordinate_system_;
};

extern template class GasPressureDependentPermeability<2>;
extern template class GasPressureDependentPermeability<3>;

}  // namespace MaterialPropertyLib
