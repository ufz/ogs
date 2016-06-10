/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_BOUNDARY_CONDITION_H_
#define PROCESS_LIB_BOUNDARY_CONDITION_H_

#include <algorithm>
#include <vector>

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h" // for GlobalIndexType

#include "DirichletBc.h"

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
public:
    UniformDirichletBoundaryCondition(GeoLib::GeoObject const* const geometry,
                                      BaseLib::ConfigTree const& config)
        : _geometry(geometry)
    {
        DBUG("Constructing UniformDirichletBoundaryCondition from config.");
        //! \ogs_file_param{boundary_condition__type}
        config.checkConfigParameter("type", "UniformDirichlet");

        //! \ogs_file_param{boundary_condition__UniformDirichlet__value}
        _value = config.getConfigParameter<double>("value");
        DBUG("Using value %g", _value);
    }

    /// Initialize Dirichlet type boundary conditions.
    /// Fills in global_ids of the particular geometry of the boundary condition
    /// and the corresponding values.
    /// The ids and the constant values are then used to construct DirichletBc
    /// object.
    void initialize(
            MeshGeoToolsLib::MeshNodeSearcher& searcher,
            NumLib::LocalToGlobalIndexMap const& dof_table,
            std::size_t component_id,
            DirichletBc<GlobalIndexType>& bc)
    {
        // Find nodes' ids on the given mesh on which this boundary condition
        // is defined.
        std::vector<std::size_t> ids = searcher.getMeshNodeIDs(*_geometry);

        // convert mesh node ids to global index for the given component
        bc.ids.reserve(bc.ids.size() + ids.size());
        bc.values.reserve(bc.values.size() + ids.size());
        for (auto& id : ids)
        {
            MeshLib::Location l(searcher.getMeshId(),
                                MeshLib::MeshItemType::Node,
                                id);
            // TODO: that might be slow, but only done once
            const auto g_idx = dof_table.getGlobalIndex(l, component_id);
            // For the DDC approach (e.g. with PETSc option), the negative
            // index of g_idx means that the entry by that index is a ghost one,
            // which should be dropped. Especially for PETSc routines MatZeroRows
            // and MatZeroRowsColumns, which are called to apply the Dirichlet BC,
            // the negative index is not accepted like other matrix or vector
            // PETSc routines. Therefore, the following if-condition is applied.
            if (g_idx >= 0)
            {
                bc.ids.emplace_back(g_idx);
                bc.values.emplace_back(_value);
            }
        }
    }

private:
    double _value;
    GeoLib::GeoObject const* const _geometry;
};


}   // namespace ProcessLib

#endif  // PROCESS_LIB_BOUNDARY_CONDITION_H_
