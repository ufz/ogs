/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "SourceTerm.h"
#include "LineSourceTermFEM.h"

namespace ProcessLib
{
class LineSourceTerm final : public SourceTerm
{
public:
    LineSourceTerm(
        unsigned const bulk_mesh_dimension,
        MeshLib::Mesh const& source_term_mesh,
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
        unsigned const integration_order, unsigned const shapefunction_order,
        ParameterLib::Parameter<double> const& source_term_parameter);

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const override;

private:
    ParameterLib::Parameter<double> const& _source_term_parameter;
    std::vector<std::unique_ptr<LineSourceTermLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace ProcessLib
