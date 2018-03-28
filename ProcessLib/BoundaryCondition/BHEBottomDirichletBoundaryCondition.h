/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BoundaryCondition.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHEAbstract.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace ProcessLib
{
class BHEBottomDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    BHEBottomDirichletBoundaryCondition(
        GlobalIndexType global_idx_T_in_bottom,
        GlobalIndexType global_idx_T_out_bottom,
        MeshLib::Mesh const& bulk_mesh,
        std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes,
        int const variable_id,
        int const component_id);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void preTimestep(const double t, const GlobalVector& x) override;

private:
    MeshLib::Mesh const& _bulk_mesh;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;

    NumLib::IndexValueVector<GlobalIndexType> _T_in_values;
};

std::unique_ptr<BHEBottomDirichletBoundaryCondition>
createBHEBottomDirichletBoundaryCondition(
    GlobalIndexType global_idx_T_in_bottom,
    GlobalIndexType global_idx_T_out_bottom, MeshLib::Mesh const& bulk_mesh,
    std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes,
    int const variable_id, int const component_id);
}  // namespace ProcessLib
