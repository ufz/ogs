/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_VARIABLE_H_
#define PROCESS_LIB_PROCESS_VARIABLE_H_

#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"
#include "ProcessLib/BoundaryCondition/BoundaryConditionConfig.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace MeshLib
{
class Mesh;
template <typename T> class PropertyVector;
}

namespace ProcessLib
{
class BoundaryConditionBuilder;

/// A named process variable. Its properties includes the mesh, and the initial
/// and boundary conditions.
class ProcessVariable
{
public:
    ProcessVariable(
        BaseLib::ConfigTree const& config, MeshLib::Mesh& mesh,
        GeoLib::GEOObjects const& geometries,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters);

    ProcessVariable(ProcessVariable&&);

    std::string const& getName() const;

    /// Returns a mesh on which the process variable is defined.
    MeshLib::Mesh const& getMesh() const;

    /// Returns the number of components of the process variable.
    int getNumberOfComponents() const { return _n_components; }

    void setBoundaryConditionBuilder(std::unique_ptr<BoundaryConditionBuilder> bc_builder)
    {
        _bc_builder = std::move(bc_builder);
    }

    std::vector<std::unique_ptr<BoundaryCondition>> createBoundaryConditions(
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
    unsigned _shapefunction_order;  ///< Order of the shapefunctions. Requires
                               /// appropriate mesh.
    Parameter<double> const& _initial_condition;

    std::vector<BoundaryConditionConfig> _bc_configs;
    std::unique_ptr<BoundaryConditionBuilder> _bc_builder;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_VARIABLE_H_
