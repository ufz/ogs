/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 4, 2020, 10:13 AM
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
 * \brief A failure index dependent permeability model \cite wangetalPerm2020
 *
 * \f[ \mathbf{k} =
 *      \mathbf{k}_0+ H(f-1) k_\text{r} \mathrm{e}^{b f}\mathbf{I}\f]
 *
 * where
 *   \f$\mathbf{k}_0\f$ is the intrinsic permeability
 *  of the undamaged material,
 * \f$H\f$ is the Heaviside step function, \f$f\f$ is the failure index,
 * \f$k_\text{r}\f$ is a reference permeability,
 * \f$b\f$ is a fitting parameter.
 * \f$k_\text{r}\f$ and \f$b\f$ can be calibrated by experimental data.
 * The failure index \f$f\f$ is calculated from  the Mohr Coulomb failure
 * criterion comparing an acting shear stress. The Mohr Coulomb failure
 * criterion \cite labuz2012mohr takes
 * the form
 *   \f[\tau(\sigma)=c-\sigma \mathrm{tan} \phi\f]
 *   with \f$\tau\f$ the shear stress, \f$c\f$ the cohesion, \f$\sigma\f$ the
 * tensile stress, and \f$\phi\f$ the internal friction angle.
 *
 *  The failure index is calculated by
 *   \f[
 *         f=\frac{\tau_m }{\cos(\phi)\tau(\sigma_m)}
 *   \f]
 *   with
 *   \f$\tau_m=(\sigma_3-\sigma_1)/2\f$
 *   and \f$\sigma_m=(\sigma_1+\sigma_3)/2\f$,
 *  where \f$\sigma_1\f$ and \f$\sigma_3\f$ are the minimum and maximum shear
 * stress, respectively.
 *
 *  Note: the conventional mechanics notations are used, which mean that tensile
 * stress is positive.
 *
 */
template <int DisplacementDim>
class PermeabilityMohrCoulombFailureIndexModel final : public Property
{
public:
    PermeabilityMohrCoulombFailureIndexModel(
        std::string name, ParameterLib::Parameter<double> const& k0,
        double const kr, double const b, double const c, double const phi,
        double const k_max,
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
    /// Intrinsic permeability for undamaged material.
    ParameterLib::Parameter<double> const& k0_;
    /// Reference permeability.
    double const kr_;
    /// Fitting parameter.
    double const b_;
    /// Cohesion.
    double const c_;
    /// Angle of internal friction.
    double const phi_;

    /// Maximum permeability.
    double const k_max_;
    ParameterLib::CoordinateSystem const* const local_coordinate_system_;
};

}  // namespace MaterialPropertyLib
