/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
        MeshLib::Mesh const& source_term_mesh,
        NumLib::LocalToGlobalIndexMap const& source_term_dof_table,
        unsigned const integration_order, unsigned const shapefunction_order,
        int const variable_id, int const component_id,
        Parameter<double> const& volumetric_source_term);

    void integrate(const double t, GlobalVector& b) const;

private:
    Parameter<double> const& _volumetric_source_term;
    std::vector<std::unique_ptr<VolumetricSourceTermLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace ProcessLib
