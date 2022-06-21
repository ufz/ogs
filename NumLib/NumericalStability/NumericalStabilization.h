/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 1, 2022, 12:14 PM
 */

#pragma once

#include <memory>
#include <vector>

namespace MeshLib
{
class Mesh;
}

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
/** It defines stabilization method for solving the advection diffusion
 * transport equation using the Galerkin finite element method.
 *
 * The convection diffusion transport equation takes the form
 *  \f[
 *     \frac{\partial u}{\partial t} - \nabla(\mathbf{K}\nabla u)
 *      + {\mathbf v}\cdot \nabla u  = Q
 *  \f]
 *  with \f$u\f$ the primary variable,  \f$\mathbf v\f$ the fluid velocity,
 *  \f$\mathbf{K}\f$ the diffusion coefficient.
 */
class NumericalStabilization
{
public:
    explicit NumericalStabilization(double const cutoff_velocity)
        : cutoff_velocity_(cutoff_velocity)
    {
    }

    virtual ~NumericalStabilization() = 0;

    virtual double getExtraDiffusionCoefficient(
        std::size_t const /*elemend_id*/,
        double const /*advection_coefficient*/,
        double const /*velocity_norm*/) const
    {
        return 0.0;
    };

protected:
    /// The cutoff velocity. The stabilization is not applied
    /// if the velocity magnitude is below the cutoff velocity.
    double const cutoff_velocity_ = 0.0;
};

/**
 * It defines an isotropic diffusion stabilization, the simplest stabilization
 * method.
 *
 *  The method adds an artificial isotropic balancing dissipation to the
 * diffusion coefficient in order to force the Péclet number to be smaller
 * than 1. The isotropic balancing dissipation is defined as
 * \f[
 *       \mathbf{K}_{\delta} = \frac{1}{2}\alpha ||\mathbf v||h \mathbf I
 *  \f]
 *  with \f$\alpha \in [0,1] \f$ the tuning parameter, \f$h\f$ the element
 *  size (e.g. the maximum edge length of element), and \f$\mathbf I\f$ the
 * identity matrix.
 */
class IsotropicDiffusionStabilization final : public NumericalStabilization
{
public:
    IsotropicDiffusionStabilization(double const cutoff_velocity,
                                    double const tuning_parameter,
                                    std::vector<double>&& element_sizes_vector);

    double getExtraDiffusionCoefficient(
        std::size_t const elemend_id,
        double const advection_coefficient,
        double const velocity_norm) const override;

private:
    /// The tuning parameter in the range [0,1].
    double const tuning_parameter_ = 0.5;

    /// Element sizes, which are represented by the maximum element edge
    /// lengths.
    std::vector<double> const element_sizes_;
};

std::unique_ptr<NumericalStabilization> createNumericalStabilization(
    MeshLib::Mesh const& mesh, BaseLib::ConfigTree const& config);
}  // namespace NumLib
