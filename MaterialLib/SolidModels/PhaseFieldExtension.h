/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct PhaseFieldExtension : public MechanicsBase<DisplacementDim>
{
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
    virtual bool calculateDegradedStress(double const t,
                                         ProcessLib::SpatialPosition const& x,
                                         KelvinVector const& eps,
                                         double& strain_energy_tensile,
                                         KelvinVector& sigma_tensile,
                                         KelvinVector& sigma_compressive,
                                         KelvinMatrix& C_tensile,
                                         KelvinMatrix& C_compressive,
                                         KelvinVector& sigma_real,
                                         double const degradation,
                                         double& elastic_energy) const = 0;

    /// Dynamic size Kelvin vector and matrix wrapper for the polymorphic
    /// constitutive relation compute function.
    bool calculateDegradedStress(
        double const t,
        ProcessLib::SpatialPosition const& x,
        Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps,
        double& strain_energy_tensile,
        Eigen::Matrix<double, Eigen::Dynamic, 1>& sigma_tensile,
        Eigen::Matrix<double, Eigen::Dynamic, 1>& sigma_compressive,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
            C_tensile,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
            C_compressive,
        Eigen::Matrix<double, Eigen::Dynamic, 1>& sigma_real,
        double const degradation,
        double& elastic_energy) const
    {
        // TODO Avoid copies of data:
        // Using MatrixBase<Derived> not possible because template functions
        // cannot be virtual. Maybe there is a workaround for this.  Using
        // Map<Matrix<double, ...>> makes the interface (for the material model
        // implementation) unnecessary difficult.
        KelvinVector const eps_{eps};
        KelvinVector sigma_tensile_{sigma_tensile};
        KelvinVector sigma_compressive_{sigma_compressive};
        KelvinMatrix C_tensile_{C_tensile};
        KelvinMatrix C_compressive_{C_compressive};
        KelvinVector sigma_real_{sigma_real};

        bool const result = calculateDegradedStress(t,
                                                    x,
                                                    eps_,
                                                    strain_energy_tensile,
                                                    sigma_tensile_,
                                                    sigma_compressive_,
                                                    C_tensile_,
                                                    C_compressive_,
                                                    sigma_real_,
                                                    degradation,
                                                    elastic_energy);

        sigma_tensile = sigma_tensile_;
        sigma_compressive = sigma_compressive_;
        C_tensile = C_tensile_;
        C_compressive = C_compressive_;
        sigma_real = sigma_real_;
        return result;
    }
};

}  // namespace Solids
}  // namespace MaterialLib
