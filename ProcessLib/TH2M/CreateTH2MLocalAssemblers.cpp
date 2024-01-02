/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTH2MLocalAssemblers.h"

#include "ProcessLib/Utils/CreateLocalAssemblersTaylorHood.h"
#include "TH2MFEM.h"

namespace ProcessLib::TH2M
{
template <int DisplacementDim>
void createLocalAssemblers(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface<DisplacementDim>>>&
        local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    bool const is_axially_symmetric,
    TH2MProcessData<DisplacementDim>& process_data)
{
    createLocalAssemblersHM<DisplacementDim, TH2MLocalAssembler>(
        mesh_elements, dof_table, local_assemblers, integration_order,
        is_axially_symmetric, process_data);
}

template void createLocalAssemblers<2>(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface<2>>>& local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    bool const is_axially_symmetric,
    TH2MProcessData<2>& process_data);

template void createLocalAssemblers<3>(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface<3>>>& local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    bool const is_axially_symmetric,
    TH2MProcessData<3>& process_data);
}  // namespace ProcessLib::TH2M
