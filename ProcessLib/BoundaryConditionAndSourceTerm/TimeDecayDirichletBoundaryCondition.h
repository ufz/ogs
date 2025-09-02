/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-07-29 11:42:08
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
 * \brief Dirichlet boundary condition with time-dependent decay.
 *
 * This boundary condition imposes the initial values of the specified primary
 * variable at the boundary nodes, scaled by a time-dependent parameter.
 * The scaling parameter should be monotonically decreasing over time,
 * representing progressive removal or reduction of support (e.g., during
 * excavation). Typically, the scaling factor starts at 1
 * (\f$ g(0, \mathbf{x}) = 1 \f$) at the initial time and decreases to 0
 * (\f$ g(t_e, \mathbf{x}) = 0 \f$) at the end of the process (\f$ t_e \f$).
 * The value of the boundary condition is given by
 * \f[ g(t, \mathbf{x}) (u_0(\mathbf{x}) - u_{\text{min}}) + u_{\text{min}} \f],
 * where \f$ u_0(\mathbf{x}) \f$ is the initial value of the primary variable
 * at the boundary node \f$ \mathbf{x} \f$, and \f$ u_{\text{min}} \f$ is the
 * user-defined lower limit of the boundary value.
 */
class TimeDecayDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    /**
     * \brief Constructs a time decay Dirichlet boundary condition.
     *
     * \param variable_id The ID of the variable to which the condition are
     *                    applied.
     * \param component_id The ID of the variable component to which the
     *                     condition are applied.
     * \param boundary_mesh The mesh representing the boundary where the nodal
     *                      forces are applied.
     * \param dof_table_bulk A unique pointer to the DOF table created from the
     *                       boundary mesh.
     * \param time_decay_parameter A parameter that defines the scaling factor
     *                             of the initial nodal values. It should be a
     *                             monotonic decrease parameter, meaning it
     *                             should decrease over time.
     * \param lower_limit The lower limit of the boundary value.
     */
    explicit TimeDecayDirichletBoundaryCondition(
        int const variable_id, int const component_id,
        MeshLib::Mesh const& boundary_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        ParameterLib::Parameter<double> const& time_decay_parameter,
        double const lower_limit);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    int const variable_id_;
    int const component_id_;

    MeshLib::Mesh const& boundary_mesh_;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap> const dof_table_boundary_;

    /** A monotonically decreasing parameter that defines the scaling factor of
     * the initial values of the specified primary variable. It is a time- and
     * position-dependent release parameter
     * \f$ g(t, \mathbf{x}) \f$
     * representing the progressive removal of support over time, e.g., \f$ g(0,
     * \mathbf{x}) = 1 \f$ and \f$ g(t_e, \mathbf{x}) = 0, \f$
     *  where \f$ t_e \f$
     * is the end time of decay, and \f$ \frac{\partial g}{\partial t} < 0
     * \f$.
     */
    ParameterLib::Parameter<double> const& time_decay_parameter_;

    /// The lower limit of the boundary value.
    double const lower_limit_;

    /// This vector is used to store the initial values at the boundary nodes
    /// before they are scaled by the time decay parameter.
    mutable std::vector<double> initial_x_values_;

    mutable bool initial_values_are_set_ = false;

    void setInitialValues(GlobalVector const& x) const;
};

}  // namespace ProcessLib
