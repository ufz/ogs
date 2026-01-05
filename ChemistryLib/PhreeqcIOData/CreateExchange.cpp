// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateExchange.h"

#include "BaseLib/ConfigTree.h"
#include "Exchange.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<ExchangeSite> createExchange(
    std::optional<BaseLib::ConfigTree> const& config, MeshLib::Mesh& mesh)
{
    if (!config)
    {
        return {};
    }

    std::vector<ExchangeSite> exchangers;
    for (auto const& site_config :
         //! \ogs_file_param{prj__chemical_system__exchangers__exchange_site}
         config->getConfigSubtreeList("exchange_site"))
    {
        //! \ogs_file_param{prj__chemical_system__exchangers__exchange_site__ion_exchanging_species}
        auto name = site_config.getConfigParameter<std::string>(
            "ion_exchanging_species");

        auto const molality = MeshLib::getOrCreateMeshProperty<double>(
            mesh, name, MeshLib::MeshItemType::IntegrationPoint, 1);

        exchangers.emplace_back(std::move(name), molality);
    }

    return exchangers;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
