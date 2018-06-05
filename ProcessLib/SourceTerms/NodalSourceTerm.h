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
#include "ProcessLib/Parameter/Parameter.h"

namespace ProcessLib
{
class NodalSourceTerm final
{
public:
    NodalSourceTerm(const NumLib::LocalToGlobalIndexMap& dof_table,
                    std::size_t const bulk_mesh_id, MeshLib::Mesh& st_mesh,
                    const int variable_id, const int component_id,
                    Parameter<double> const& parameter);

    void integrateNodalSourceTerm(const double t, GlobalVector& b) const;

private:
    NumLib::LocalToGlobalIndexMap const& _dof_table;
    std::size_t const _bulk_mesh_id;
    MeshLib::Mesh& _st_mesh;
    int const _variable_id;
    int const _component_id;
    Parameter<double> const& _parameter;
};

}  // namespace ProcessLib
