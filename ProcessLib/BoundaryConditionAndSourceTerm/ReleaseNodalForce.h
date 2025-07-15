/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-07-04 11:34:09
 */

#pragma once

#include <memory>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/BoundaryCondition.h"

namespace ParameterLib
{
template <typename T>
struct Parameter;
}  // namespace ParameterLib

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{

/**
 * \brief Boundary condition for simulating excavation using the release nodal
 * force approach.
 *
 * This class implements a boundary condition that applies a time-dependent,
 * released nodal force to nodes on an exposed surface for excavation
 * simulations.
 *
 * The initial state assumes a non-equilibrium stress \f$ \sigma_0 \f$. In the
 * finite element method, the nodal force is given by:
 * \f[
 *   \mathbf{b} = \int \left( \text{B}^T (\mathbf{\sigma} - \mathbf{\sigma}_0) +
 *   (\mathbf{f} - \mathbf{f}_0) \right) N \, \mathrm{d}\Omega + \int_{\Gamma_q}
 * (\boldsymbol{\sigma}-\boldsymbol{\sigma}_0)\cdot \mathbf n \mathrm{d}\Gamma
 * \f]
 * where:
 *   - \f$ \text{B} \f$ is the strain-displacement matrix,
 *   - \f$ \mathbf{\sigma} \f$ is the current total stress,
 *   - \f$ \mathbf{\sigma}_0 \f$ is the initial total stress,
 *   - \f$ \mathbf{f} \f$ is the current body force,
 *   - \f$ \mathbf{f}_0 \f$ is the initial body force,
 *   - \f$\int_{\Gamma_q}\f$ is the boundary where the traction condition is
 * applied, and \f$ \mathbf n\f$ is the outward normal of the boundary,
 *   - \f$ N \f$ is the shape function,
 *   - \f$ \Omega \f$ is the domain of integration.
 *
 * After excavation, the stress and body force inside the excavated domain
 * vanish, leaving non-zero nodal forces at the exposed surface nodes. These are
 * computed as:
 * \f[
 *   \mathbf{b}_0 = -\int \left( \text{B}^T \mathbf{\sigma}_0 + \mathbf{f}_0 N
 * \right)  \, \mathrm{d}\Omega - \int_{\Gamma_q} \boldsymbol{\sigma}_0 \cdot
 * \mathbf n \mathrm{d}\Gamma
 * \f]
 * where \f$\Omega\f$ is the remaining domain.
 *
 * The elements of \f$ \mathbf{b}_0 \f$ corresponding to the exposed surface
 * nodes define the released nodal force vector:
 * \f[
 *   \mathbf{f}_\text{r} := (\mathbf{b}_0)_i, \quad i \in \text{exposed surface
 * nodes}.
 * \f]
 *  \f$\Omega\f$ can be the excavated domain, which leads to the negative \f$
 * \mathbf{f}_\text{r}\f$.
 *
 * To simulate excavations under an assumption of gradual release of these
 * forces, the boundary condition applies the released nodal force vector to the
 * global right-hand side (RHS) vector \c b, scaled by a time- and
 * position-dependent release parameter \f$ g(t, \mathbf{x}) \f$:
 * \f[
 *   \mathbf{b} = \mathbf{b} + \mathbf{f}_\text{r} \cdot g(t, \mathbf{x})
 * \f]
 * The release parameter should be a monotonically decreasing function,
 * representing the progressive removal of support over time, e.g., \f$ g(0,
 * \mathbf{x}) = 1 \f$ and \f$ g(t_e, \mathbf{x}) = 0 \f$, where \f$ t_e \f$ is
 * the end time of excavation, and \f$ \frac{\partial g}{\partial t} < 0 \f$.
 *
 * This boundary condition is particularly useful for modeling staged
 * excavations or similar processes where loads are released in a controlled
 * manner over time.
 *
 * \note Setting `compensate_non_equilibrium_initial_residuum` to true in the
 * process variable configuration is required when using this boundary condition
 * to ensure that the initial non-equilibrium stress state is properly accounted
 * for.
 */
class ReleaseNodalForce final : public BoundaryCondition
{
public:
    /**
     * \brief Constructs a released nodal force boundary condition.
     *
     * \param variable_id The ID of the variable to which the nodal forces are
     *                    applied.
     * \param boundary_mesh The mesh representing the boundary where the nodal
     *                      forces are applied.
     * \param dof_table A unique pointer to the DOF table created from the
     *                  boundary mesh.
     * \param time_decay_parameter A parameter that defines the scaling factor
     *                             of the nodal force. It should be a monotonic
     *                             decrease parameter, meaning it should
     *                             decrease over time.
     */
    explicit ReleaseNodalForce(
        int const variable_id,
        MeshLib::Mesh const& boundary_mesh,
        std::unique_ptr<NumLib::LocalToGlobalIndexMap>& dof_table,
        ParameterLib::Parameter<double> const& time_decay_parameter);

    void set(GlobalVector const* r_neq);

    /**
     * \brief Applies the released nodal force boundary condition.
     * This method scales the nodal forces by the release parameter and
     * applies them to the global RHS vector \c b.
     * \param t The current time.
     * \param x The global solution vector (not used in this case).
     * \param process_id The ID of the process.
     * \param K The global stiffness matrix (not used in this case).
     * \param b The global RHS vector to which the nodal forces are applied.
     * \param Jac The global Jacobian matrix (not used in this case).
     */
    void applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                        int const process_id, GlobalMatrix* K, GlobalVector& b,
                        GlobalMatrix* Jac) override;

private:
    int const variable_id_;

    MeshLib::Mesh const& boundary_mesh_;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap> const dof_table_;

    /** A monotonically decreasing parameter that defines the scaling factor of
     * the nodal force. It is a time- and position-dependent release parameter
     * \f$ g(t, \mathbf{x}) \f$
     * representing the progressive removal of support over time, e.g., \f$ g(0,
     * \mathbf{x}) = 1 \f$ and \f$ g(t_e, \mathbf{x}) = 0 \f$, where \f$ t_e \f$
     * is the end time of excavation, and \f$ \frac{\partial g}{\partial t} < 0
     * \f$.
     */
    ParameterLib::Parameter<double> const& time_decay_parameter_;

    /// Initial nodal force values at the boundary nodes.
    /// This vector is used to store the initial nodal forces before they are
    /// scaled by the release parameter.
    std::vector<double> initial_release_nodal_force_;

    std::vector<GlobalIndexType> global_indices_;
    std::vector<MeshLib::Node const*> boundary_nodes_;
};

}  // namespace ProcessLib
