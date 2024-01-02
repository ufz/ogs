/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Element.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/IntegrationMethodRegistry.h"

namespace ProcessLib::TH2M
{
template <int DisplacementDim>
struct LocalAssemblerInterface;

template <int DisplacementDim>
struct TH2MProcessData;

template <int DisplacementDim>
void createLocalAssemblers(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface<DisplacementDim>>>&
        local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    bool const is_axially_symmetric,
    TH2MProcessData<DisplacementDim>& process_data);
}  // namespace ProcessLib::TH2M
