/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/NumericsConfig.h"

namespace GeoLib
{
class GeoObject;
}

namespace MeshLib
{
class Element;
class Mesh;
}

namespace MeshGeoToolsLib
{
class MeshNodeSearcher;
class BoundaryElementsSearcher;
}

namespace NumLib
{
class LocalToGlobalIndexMap;
template <typename>
struct IndexValueVector;
}

namespace ProcessLib
{
struct BoundaryConditionConfig;
struct ParameterBase;

class BoundaryCondition
{
public:
    //! Applies natural BCs (i.e. non-Dirichlet BCs) to the stiffness matrix
    //! \c K and the vector \c b.
    virtual void applyNaturalBC(const double /*t*/, GlobalVector const& /*x*/,
                                GlobalMatrix& /*K*/, GlobalVector& /*b*/)
    {
        // By default it is assumed that the BC is not a natural BC. Therefore
        // there is nothing to do here.
    }

    //! Writes the values of essential BCs to \c bc_values.
    virtual void getEssentialBCValues(
        const double /*t*/,
        NumLib::IndexValueVector<GlobalIndexType>& /*bc_values*/) const
    {
        // By default it is assumed that the BC is not an essential BC.
        // Therefore there is nothing to do here.
    }

    virtual void preTimestep(const double /*t*/) {}

    virtual ~BoundaryCondition() = default;
};

class BoundaryConditionBuilder
{
public:
    virtual ~BoundaryConditionBuilder() {}

    virtual std::unique_ptr<BoundaryCondition> createBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
        const int variable_id, const unsigned integration_order,
        const unsigned shapefunction_order,
        const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

protected:
    virtual std::unique_ptr<BoundaryCondition> createDirichletBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const MeshLib::Mesh& mesh, const int variable_id,
        const unsigned integration_order,
        const unsigned shapefunction_order,
        const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>&
            parameters,
        MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher,
        MeshGeoToolsLib::BoundaryElementsSearcher& boundary_element_searcher);

    virtual std::unique_ptr<BoundaryCondition> createNeumannBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const MeshLib::Mesh& mesh, const int variable_id,
        const unsigned integration_order,
        const unsigned shapefunction_order,
        const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>&
            parameters,
        MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher,
        MeshGeoToolsLib::BoundaryElementsSearcher& boundary_element_searcher);

    virtual std::unique_ptr<BoundaryCondition> createRobinBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const MeshLib::Mesh& mesh, const int variable_id,
        const unsigned integration_order,
        const unsigned shapefunction_order,
        const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>&
            parameters,
        MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher,
        MeshGeoToolsLib::BoundaryElementsSearcher& boundary_element_searcher);

    static std::vector<MeshLib::Element*> getClonedElements(
        MeshGeoToolsLib::BoundaryElementsSearcher& boundary_element_searcher,
        GeoLib::GeoObject const& geometry);
};

}  // ProcessLib
