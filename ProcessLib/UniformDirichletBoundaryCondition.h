/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_BOUNDARY_CONDITION_H_
#define PROCESS_LIB_BOUNDARY_CONDITION_H_

#include <algorithm>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "AssemblerLib/MeshComponentMap.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

#include "MeshLib/Mesh.h"

namespace GeoLib
{
    class GeoObject;
}

namespace ProcessLib
{

/// The UniformDirichletBoundaryCondition class describes a constant in space
/// and time Dirichlet boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class UniformDirichletBoundaryCondition
{
    using ConfigTree = boost::property_tree::ptree;
public:
    UniformDirichletBoundaryCondition(GeoLib::GeoObject const* const geometry,
            ConfigTree const& config)
        : _geometry(geometry)
    {
        DBUG("Constructing UniformDirichletBoundaryCondition from config.");

        _value = config.get<double>("value", 0);
        DBUG("Using value %g", _value);
    }

    /// Initialize Dirichlet type boundary conditions.
    /// Fills in global_ids of the particular geometry of the boundary condition
    /// and the corresponding values.
    /// The ids are appended to the global_ids and the values are filled with
    /// the constant _value.
    void initialize(
            MeshGeoToolsLib::MeshNodeSearcher& searcher,
            AssemblerLib::MeshComponentMap const& mcmap,
            std::size_t component_id,
            std::vector<std::size_t>& global_ids,
            std::vector<double>& values)
    {
        // Find nodes' ids on the given mesh on which this boundary condition
        // is defined.
        std::vector<std::size_t> ids = searcher.getMeshNodeIDs(*_geometry);

        // convert mesh node ids to global index for the given component
        for (auto& id : ids)
        {
            id = mcmap.getGlobalIndex(
                     MeshLib::Location(searcher.getMesh().getID(),
                                       MeshLib::MeshItemType::Node,
                                       id),
                     component_id);
            // id = nodal_dof_count * id + nodal_dof_index; // order by mesh node
            // id = mesh_node_count * nodal_dof_index + id; // order by variable
        }

        // Append node ids.
        global_ids.reserve(global_ids.size() + ids.size());
        std::copy(ids.cbegin(), ids.cend(), std::back_inserter(global_ids));

        // Fill values.
        values.reserve(values.size() + ids.size());
        std::fill_n(std::back_inserter(values), ids.size(), _value);
    }

private:
    double _value;
    GeoLib::GeoObject const* const _geometry;
};


}   // namespace ProcessLib

#endif  // PROCESS_LIB_BOUNDARY_CONDITION_H_
