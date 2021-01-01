/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
        std::size_t const source_term_mesh_id, MeshLib::Mesh const& st_mesh,
        const int variable_id, const int component_id,
        ParameterLib::Parameter<double> const& parameter);

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const override;

private:
    std::size_t const _source_term_mesh_id;
    MeshLib::Mesh const& _st_mesh;
    int const _variable_id;
    int const _component_id;
    ParameterLib::Parameter<double> const& _parameter;
};

}  // namespace ProcessLib
