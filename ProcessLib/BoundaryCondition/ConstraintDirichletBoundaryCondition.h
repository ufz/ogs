/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "BoundaryCondition.h"

#include "ConstraintDirichletBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
/// The ConstraintDirichletBoundaryCondition class describes a Dirichlet-type
/// boundary condition that is constant in space and time where the domain can
/// shrink and grow within the simulation. The expected parameter in the passed
/// configuration is "value" which, when not present defaults to zero.
class ConstraintDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    /// @param parameter Used for setting the values for the boundary condition.
    /// @param dof_table_bulk The bulk local to global index map is used to
    /// derive the local to global index map for the boundary.
    /// @param variable_id The variable id is needed to determine the global
    /// index.
    /// @param component_id The component id is needed to determine the global
    /// index.
    /// @param bc_mesh Lower dimensional mesh the boundary condition is defined
    /// on. The bc_mesh must have the two PropertyVector objects
    /// 'bulk_element_ids' and 'bulk_node_ids' containing the corresponding
    /// information.
    /// @param integration_order Order the order of integration used to compute
    /// the constraint.
    /// @param bulk_mesh The FE mesh for the simulation.
    /// @param constraint_threshold The threshold value used for the switch
    /// off/on decision.
    /// @param lower Boolean value used for the calculation of the constraint
    /// criterion, i.e., if lower is set to true the criterion 'calculated_value
    /// < constraint_threshold' is evaluated to switch on/off the boundary
    /// condition, else 'calculated_value > constraint_threshold' is evaluated.
    /// @param getFlux The function used for the flux calculation.
    /// @note The function has to be stored by value, else the process value is
    /// not captured properly.
    ConstraintDirichletBoundaryCondition(
        Parameter<double> const& parameter,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        MeshLib::Mesh const& bc_mesh, unsigned const integration_order,
        MeshLib::Mesh const& bulk_mesh, double const constraint_threshold,
        bool const lower,
        std::function<Eigen::Vector3d(std::size_t const,
                                      MathLib::Point3d const&, double const,
                                      GlobalVector const&)>
            getFlux);

    void preTimestep(double t, GlobalVector const& x) override;

    void getEssentialBCValues(
        const double t,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    Parameter<double> const& _parameter;

    /// Local dof table, a subset of the global one restricted to the
    /// participating number of elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    int const _variable_id;
    int const _component_id;

    /// Vector of (lower-dimensional) boundary elements on which the boundary
    /// condition is defined.
    MeshLib::Mesh const& _bc_mesh;

    /// Integration order for integration over the lower-dimensional elements
    unsigned const _integration_order;

    /// The first item of the pair is the element id in the bulk mesh, the
    /// second item is the face id of the bulk element that is part of the
    /// boundary
    std::vector<std::pair<std::size_t, unsigned>> _bulk_ids;

    /// Stores the results of the flux computations per boundary element.
    std::vector<double> _flux_values;

    /// Local assemblers for each boundary element.
    std::vector<std::unique_ptr<
        ConstraintDirichletBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;

    /// The threshold value used to the switch off/on the Dirichlet-type
    /// boundary condition.
    double const _constraint_threshold;

    /// The boolean value lower is used for the calculation of the constraint
    /// criterion, i.e., if lower is set to true the criterion 'calculated_value
    /// < constraint_threshold' is evaluated to switch on/off the boundary
    /// condition, else 'calculated_value > constraint_threshold' is evaluated.
    bool const _lower;

    /// The mesh _bulk_mesh is the discretized domain the process(es) are
    /// defined on. It is needed to get values for the constraint calculation.
    MeshLib::Mesh const& _bulk_mesh;

    /// The function _getFlux calculates the flux through the boundary element.
    std::function<Eigen::Vector3d(
            std::size_t const, MathLib::Point3d const&, double const,
            GlobalVector const&)> _getFlux;
};

/// The function parses the config tree and creates a
/// ConstraintDirichletBoundaryCondition.
std::unique_ptr<ConstraintDirichletBoundaryCondition>
createConstraintDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    unsigned const integration_order, int const component_id,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    Process const& constraining_process);
}  // namespace ProcessLib
