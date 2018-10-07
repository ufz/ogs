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
class BHEInflowDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    BHEInflowDirichletBoundaryCondition(
        std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
        MeshLib::Mesh const& bc_mesh,
        std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
        int const variable_id,
        int const component_id,
        std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
            pt_bhe);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void preTimestep(const double t, const GlobalVector& x) override;

private:
    // TODO (haibing) re-organize as the bottom BC data structure
    MeshLib::Mesh const& _bc_mesh;

    /// Stores the results of the outflow temperatures per boundary node.
    std::vector<double> _T_out_values;
    std::vector<GlobalIndexType> _T_out_indices;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;

    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
        _pt_bhe;
};

std::unique_ptr<BHEInflowDirichletBoundaryCondition>
createBHEInflowDirichletBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes,
    int const variable_id, int const component_id,
    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
        pt_bhe);
}  // namespace ProcessLib
