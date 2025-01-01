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

#include "SourceTerm.h"

namespace ProcessLib
{
template <int GlobalDim>
class AnchorTerm final : public SourceTerm
{
public:
    explicit AnchorTerm(
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
        std::size_t const source_term_mesh_id, MeshLib::Mesh const& st_mesh,
        const int variable_id,
        ParameterLib::Parameter<double> const& parameter);

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const override;

private:
    std::size_t const source_term_mesh_id_;
    MeshLib::Mesh const& st_mesh_;
    int const variable_id_;
    ParameterLib::Parameter<double> const& parameter_;
};

extern template class AnchorTerm<2>;
extern template class AnchorTerm<3>;
}  // namespace ProcessLib
