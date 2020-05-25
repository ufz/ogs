/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

#include "ParameterLib/Parameter.h"
#include "ProcessLib/BoundaryCondition/BoundaryConditionConfig.h"
#include "ProcessLib/SourceTerms/SourceTermConfig.h"

// DeactivatedSubdomain cannot be forwardly declared because that
// std::unique_ptr<DeactivatedSubdomain> type member requires its full
// definition (see https://stackoverflow.com/a/6089065).
#include "ProcessLib/DeactivatedSubdomain.h"

namespace MeshLib
{
class Mesh;
template <typename T>
class PropertyVector;
}
namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
class SourceTerm;
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
        BaseLib::ConfigTree const& config, MeshLib::Mesh& mesh,
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters);

    ProcessVariable(ProcessVariable&& other);

    std::string const& getName() const;

    /// Returns a mesh on which the process variable is defined.
    MeshLib::Mesh const& getMesh() const;

    std::vector<std::unique_ptr<DeactivatedSubdomain const>> const&
    getDeactivatedSubdomains() const
    {
        return deactivated_subdomains_;
    }

    void updateDeactivatedSubdomains(double const time);

    std::vector<std::size_t>& getActiveElementIDs() const
    {
        return ids_of_active_elements_;
    }

    /// Returns the number of components of the process variable.
    int getNumberOfComponents() const { return n_components_; }
    std::vector<std::unique_ptr<BoundaryCondition>> createBoundaryConditions(
        const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
        unsigned const integration_order,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        Process const& process);

    std::vector<std::unique_ptr<SourceTerm>> createSourceTerms(
        const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
        unsigned const integration_order,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters);

    ParameterLib::Parameter<double> const& getInitialCondition() const
    {
        return initial_condition_;
    }

    unsigned getShapeFunctionOrder() const { return shapefunction_order_; }
private:
    std::string const name_;
    MeshLib::Mesh& mesh_;
    const int n_components_;
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
    unsigned shapefunction_order_;

    std::vector<std::unique_ptr<DeactivatedSubdomain const>>
        deactivated_subdomains_;

    /// IDs of the active elements. It is initialized only if there are
    /// deactivated subdomains.
    mutable std::vector<std::size_t> ids_of_active_elements_;

    void createBoundaryConditionsForDeactivatedSubDomains(
        const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        std::vector<std::unique_ptr<BoundaryCondition>>& bcs);

    ParameterLib::Parameter<double> const& initial_condition_;

    std::vector<BoundaryConditionConfig> bc_configs_;
    std::vector<SourceTermConfig> source_term_configs_;
};

}  // namespace ProcessLib
