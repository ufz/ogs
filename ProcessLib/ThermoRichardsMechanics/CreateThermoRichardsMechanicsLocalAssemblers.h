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

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim, typename ConstitutiveTraits>
struct LocalAssemblerInterface;

template <int DisplacementDim, typename ConstitutiveTraits>
struct ThermoRichardsMechanicsProcessData;

template <int DisplacementDim, typename ConstitutiveTraits>
void createLocalAssemblers(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<
        LocalAssemblerInterface<DisplacementDim, ConstitutiveTraits>>>&
        local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    bool const is_axially_symmetric,
    ThermoRichardsMechanicsProcessData<DisplacementDim, ConstitutiveTraits>&
        process_data);
}  // namespace ProcessLib::ThermoRichardsMechanics
