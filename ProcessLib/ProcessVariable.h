/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "BaseLib/ConfigTree-fwd.h"
#include "MaterialLib/MPL/Medium.h"

namespace MathLib
{
class PiecewiseLinearInterpolation;
}  // namespace MathLib

namespace MeshLib
{
class Mesh;

template <typename PROP_VAL_TYPE>
class PropertyVector;
}  // namespace MeshLib

namespace NumLib
{
class LocalToGlobalIndexMap;
}  // namespace NumLib

namespace ParameterLib
{
struct ParameterBase;

template <typename T>
struct Parameter;
}  // namespace ParameterLib

namespace ProcessLib
{
class SourceTermBase;
class BoundaryCondition;
class Process;
struct DeactivatedSubdomain;
struct BoundaryConditionConfig;
struct SourceTermConfig;
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
            parameters,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves);

    ProcessVariable(ProcessVariable&&);

    std::string const& getName() const;

    /// Returns a mesh on which the process variable is defined.
    MeshLib::Mesh const& getMesh() const;

    std::vector<DeactivatedSubdomain> const& getDeactivatedSubdomains() const
    {
        return _deactivated_subdomains;
    }

    void updateDeactivatedSubdomains(double const time);

    std::vector<std::size_t> const& getActiveElementIDs() const
    {
        return _ids_of_active_elements;
    }

    /// Returns the number of components of the process variable.
    int getNumberOfGlobalComponents() const { return _n_components; }

    std::vector<std::unique_ptr<BoundaryCondition>> createBoundaryConditions(
        const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
        unsigned const integration_order,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        Process const& process,
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            all_process_variables_for_this_process,
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    std::vector<std::unique_ptr<SourceTermBase>> createSourceTerms(
        const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
        unsigned const integration_order,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        std::vector<std::reference_wrapper<ProcessVariable>> const&
            all_process_variables_for_this_process,
        const MeshLib::Mesh& bulk_mesh);

    ParameterLib::Parameter<double> const& getInitialCondition() const
    {
        return _initial_condition;
    }

    unsigned getShapeFunctionOrder() const { return _shapefunction_order; }

    bool compensateNonEquilibriumInitialResiduum() const
    {
        return _compensate_non_equilibrium_initial_residuum;
    }

    ~ProcessVariable();

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

    std::vector<DeactivatedSubdomain> _deactivated_subdomains;

    /// IDs of the active elements. It is initialized only if there are
    /// deactivated subdomains.
    mutable std::vector<std::size_t> _ids_of_active_elements;
    MeshLib::PropertyVector<unsigned char>* _is_active = nullptr;

    void createBoundaryConditionsForDeactivatedSubDomains(
        const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        std::vector<std::unique_ptr<BoundaryCondition>>& bcs);

    ParameterLib::Parameter<double> const& _initial_condition;

    std::vector<BoundaryConditionConfig> _bc_configs;
    std::vector<SourceTermConfig> _source_term_configs;

    const bool _compensate_non_equilibrium_initial_residuum;
};

}  // namespace ProcessLib
