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

namespace ChemistryLib
{
namespace PhreeqcIOData
{
AqueousSolution createAqueousSolution(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__chemical_system__solution__temperature}
    auto const temperature = config.getConfigParameter<double>("temperature");

    //! \ogs_file_param{prj__chemical_system__solution__pressure}
    auto const pressure = config.getConfigParameter<double>("pressure");

    //! \ogs_file_param{prj__chemical_system__solution__pe}
    auto const pe = config.getConfigParameter<double>("pe");

    auto components = createSolutionComponents(config);

    auto charge_balance = createChargeBalance(config);

    return {temperature, pressure, pe, std::move(components), charge_balance};
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
