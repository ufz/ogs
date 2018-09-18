/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "SourceTerm.h"

namespace ProcessLib
{
class NodalSourceTerm final : public SourceTerm
{
public:
    explicit NodalSourceTerm(
        const NumLib::LocalToGlobalIndexMap& source_term_dof_table,
        std::size_t const bulk_mesh_id, MeshLib::Mesh const& st_mesh,
        const int variable_id, const int component_id,
        Parameter<double> const& parameter);

    void integrate(const double t, GlobalVector& b) const override;

private:
    std::size_t const _bulk_mesh_id;
    MeshLib::Mesh const& _st_mesh;
    int const _variable_id;
    int const _component_id;
    Parameter<double> const& _parameter;
};

}  // namespace ProcessLib
