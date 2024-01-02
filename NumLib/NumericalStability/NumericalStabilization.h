/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 1, 2022, 12:14 PM
 */

#pragma once

#include <Eigen/Core>
#include <variant>
#include <vector>

#include "BaseLib/Error.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"

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
struct NoStabilization final
{
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
class IsotropicDiffusionStabilization final
{
public:
    IsotropicDiffusionStabilization(double const cutoff_velocity,
                                    double const tuning_parameter,
                                    std::vector<double>&& element_sizes);

    double computeArtificialDiffusion(std::size_t const element_id,
                                      double const velocity_norm) const;

private:
    /// The cutoff velocity. The stabilization is not applied
    /// if the velocity magnitude is below the cutoff velocity.
    double const cutoff_velocity_;

    /// The tuning parameter in the range [0,1].
    double const tuning_parameter_ = 0.5;

    /// Element sizes, which are represented by the maximum element edge
    /// lengths.
    std::vector<double> const element_sizes_;
};

/**
 *
 *  For the full upwind scheme, we consider the general advection term of the
 *  advection  diffusion transport equation:
 *  \f[
 *          \nabla \cdot ( u \mathbf{v}),
 *   \f]
 *     where \f$u\f$ can be temperature \f$ T\f$ or mass component
 *  concentration \f$C\f$, and \f$\mathbf{v}\f$ is the fluid velocity.
 *  This means \f$\nabla \cdot \mathbf{v}\f$ may not be zero, or physically
 *  the fluid flow may be compressible.
 *
 *    The discretized weak form of that advection term takes the form
 *    \f[
 *        \int_{\Omega_e}  \nabla \cdot ( u \mathbf{v})
 *                \phi_i \mathrm{d} \Omega
 *              =  \int_{\Omega_e}  \nabla \cdot ( u \mathbf{v} \phi_i)
 *                 \mathrm{d} \Omega - \int_{\Omega_e} \nabla \phi_i
 *                  \cdot ( u \mathbf{v} )  \mathrm{d} \Omega
 *    \f]
 *   with \f$\phi_i\f$ the test function.
 *    The first term on the right hand side can be converted to boundary
 *  integration according to the Green's theorem as
 *          \f[
 *               \int_{\Omega_e}  \nabla \cdot ( u \mathbf{v} \phi_i)
 *                 \mathrm{d} \Omega
 *               = \int_{\partial\Omega_e}  ( u \mathbf{v} \phi_i) \cdot
 *                    \mathbf{n} \mathrm{d} \Gamma
 *          \f]
 *        with \f$\mathbf{n}\f$ the normal to the boundary surface. That
 *        boundary integration is a part of Neumann boundary condition.
 *     Therefore what left for the element integration is
 *       \f[
 *            -\int_{\Omega_e} \nabla \phi_i
 *                  \cdot ( u \mathbf{v} )  \mathrm{d} \Omega,
 *        \f]
 * which is denoted as \f$R_i\f$ hereafter.
 *
 *  Based on the scheme introduced by Dalen \cite dalen1979, the full upwind
 scheme
 *  evaluates the following flux related quantity for each node
 * \f[
 *       q_i = -\int_{\Omega_e} \nabla \phi_i
 *                  \cdot  \mathbf{v}  \mathrm{d} \Omega
 * \f]
 *  to determine whether it is an upwind node. If \f$q_i>0\f$, node \f$ i\f$ is
 * at upwind position.
 *
 *  Let
 *     \f[
 *          {q}_{up}=\sum_{q_i \geq 0} q_i u_i
 *     \f]
 *    be the total flux at the upwind nodes, and
 *     \f[
 *         {q}_{down}=-\sum_{q_i < 0} q_i
 *     \f]
 *     be the total flux related quantity at the down nodes, we can approximate
 *  the discretized weak form of the advection term as:
 *  \f{eqnarray*}{
 *         R_i
 *         & \approx&
 *           \begin{cases}
 *             q_i\,u_i,\forall q_i \geq 0\\
 *             q_i \frac{{q}_{up}}{{q}_{down}},\forall q_i < 0
 *           \end{cases}\\
 *         & = &
 *           \begin{cases}
 *             q_i\,u_i,\forall q_i \geq 0\\
 *             \frac{1}{{q}_{down}} q_i {\sum_{q_j \geq 0} (q_j u_j)},\forall
 *             q_i < 0
 *           \end{cases} =\tilde{R}_i.\\
 *        \f}
 *  The above approximation defines the full upwind scheme of the advection term
 *  of the advection diffusion transport equation for the FEM analysis.
 *  Obviously, we see that
 * \f{eqnarray*}{
 *   \sum_i   \tilde{R}_i &=& \sum_{q_i \geq 0} q_i u_i + \sum_{q_i < 0}
 *  \frac{1}{{q}_{down}} q_i {\sum_{q_j \geq 0} (q_j u_j)}\\
 *     &=& \sum_{q_i \geq 0} q_i u_i +
 *  \frac{\sum_{q_i < 0} q_i}{{q}_{down}}  {\sum_{q_j \geq 0} (q_j u_j)}
 *   = q_{up}-q_{up} = 0,
 * \f}
 * which means the nodal mass balance of element is guaranteed.
 *
 */
class FullUpwind final
{
public:
    explicit FullUpwind(double const cutoff_velocity)
        : cutoff_velocity_(cutoff_velocity)
    {
    }

    double getCutoffVelocity() const { return cutoff_velocity_; }

private:
    /// The cutoff velocity. The stabilization is not applied
    /// if the velocity magnitude is below the cutoff velocity.
    double const cutoff_velocity_;
};

class FluxCorrectedTransport final
{
public:
    explicit FluxCorrectedTransport()
    {
#ifdef USE_PETSC
        OGS_FATAL(
            "FluxCorrectedTransport scheme is not implemented to work with MPI "
            "parallelization.");
#endif
    }

    void calculateFluxCorrectedTransport(
        const double t, const double dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        const MathLib::MatrixSpecifications& matrix_specification,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b);
};

using NumericalStabilization =
    std::variant<NoStabilization, IsotropicDiffusionStabilization, FullUpwind,
                 FluxCorrectedTransport>;
}  // namespace NumLib
