/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateAqueousSolution.h"
#include "AqueousSolution.h"
#include "BaseLib/ConfigTree.h"
#include "ChemistryLib/Common/CreateChargeBalance.h"
#include "CreateSolutionComponent.h"
#include "MeshLib/Mesh.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::unique_ptr<AqueousSolution> createAqueousSolution(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh,
    MeshLib::PropertyVector<std::size_t> const& chemical_system_map)
{
    //! \ogs_file_param{prj__chemical_system__solution__temperature}
    auto const temperature = config.getConfigParameter<double>("temperature");

    //! \ogs_file_param{prj__chemical_system__solution__pressure}
    auto const pressure = config.getConfigParameter<double>("pressure");

    //! \ogs_file_param{prj__chemical_system__solution__pe}
    auto const pe0 = config.getConfigParameter<double>("pe");

    auto pe = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "pe", MeshLib::MeshItemType::Node, 1);

    std::fill(std::begin(*pe),
              std::end(*pe),
              std::numeric_limits<double>::quiet_NaN());

    std::for_each(
        chemical_system_map.begin(),
        chemical_system_map.end(),
        [&pe, pe0](auto const& global_id) { (*pe)[global_id] = pe0; });

    auto components =
        createSolutionComponents(config, mesh.getNumberOfBaseNodes());

    auto charge_balance = createChargeBalance(config);

    return std::make_unique<AqueousSolution>(temperature,
                                             pressure,
                                             pe,
                                             std::move(components),
                                             charge_balance,
                                             mesh.getNumberOfBaseNodes());
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
