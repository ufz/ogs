// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "SourceTerm.h"
#include "VolumetricSourceTermFEM.h"

namespace ProcessLib
{
class VolumetricSourceTerm final : public SourceTerm
{
public:
    VolumetricSourceTerm(
        unsigned const bulk_mesh_dimension,
        MeshLib::Mesh const& source_term_mesh,
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
        unsigned const integration_order, unsigned const shapefunction_order,
        ParameterLib::Parameter<double> const& source_term_parameter);

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const override;

private:
    ParameterLib::Parameter<double> const& _source_term_parameter;
    std::vector<std::unique_ptr<VolumetricSourceTermLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace ProcessLib
