/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/BoundaryCondition/BoundaryConditionConfig.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/SourceTerms/SourceTermConfig.h"

namespace MeshLib
{
class Mesh;
template <typename T> class PropertyVector;
}
namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
class NodalSourceTerm;
class BoundaryCondition;
class Process;
}  // namespace ProcessLib

namespace ProcessLib
{
/// A named process variable. Its properties includes the mesh, and the initial
/// and boundary conditions as well as the source terms.
class ProcessVariable
{
public:
    ProcessVariable(
        BaseLib::ConfigTree const& config,
        MeshLib::Mesh& mesh,
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters);

    ProcessVariable(ProcessVariable&&);

    std::string const& getName() const;

    /// Returns a mesh on which the process variable is defined.
    MeshLib::Mesh const& getMesh() const;

    /// Returns the number of components of the process variable.
    int getNumberOfComponents() const { return _n_components; }

    std::vector<std::unique_ptr<BoundaryCondition>> createBoundaryConditions(
        const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
        unsigned const integration_order,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        Process const& process);

    std::vector<std::unique_ptr<NodalSourceTerm>> createSourceTerms(
        const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
        unsigned const integration_order,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters);

    Parameter<double> const& getInitialCondition() const
    {
        return _initial_condition;
    }

    // Get or create a property vector for results.
    // The returned mesh property size is number of mesh nodes times number of
    // components.
    MeshLib::PropertyVector<double>& getOrCreateMeshProperty();

    unsigned getShapeFunctionOrder() const { return _shapefunction_order; }

private:
    std::string const _name;
    MeshLib::Mesh& _mesh;
    const int _n_components;
    /// The polynomial order of the process variable's shape functions.
    ///
    /// Requires an appropriate mesh.
    ///
    /// The order of the shape functions can not be higher than the maximum
    /// available order for the underlying geometric elements. For example the
    /// second order shape functions for a hexahedron are only possible if the
    /// geometric element is at least a 20-node hexahedron element
    /// (MeshLib::TemplateElement<MeshLib::HexRule20>), whereas linear shape
    /// functions are also available on the 8-node hexahedron
    /// (MeshLib::TemplateElement<MeshLib::HexRule8>).
    ///
    /// \sa MeshLib::CellRule MeshLib::FaceRule MeshLib::EdgeRule.
    unsigned _shapefunction_order;
    Parameter<double> const& _initial_condition;

    std::vector<BoundaryConditionConfig> _bc_configs;
    std::vector<SourceTermConfig> _source_term_configs;
};

}  // namespace ProcessLib
