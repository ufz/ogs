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

    ParameterLib::Parameter<double> const& time_decay_parameter_;

    /// Initial nodal force values at the boundary nodes.
    /// This vector is used to store the initial nodal forces before they are
    /// scaled by the release parameter.
    std::vector<double> initial_release_nodal_force_;

    std::vector<GlobalIndexType> global_indices_;
    std::vector<MeshLib::Node const*> boundary_nodes_;
};

}  // namespace ProcessLib
