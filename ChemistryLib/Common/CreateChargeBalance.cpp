/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateChargeBalance.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "ChargeBalance.h"

namespace ChemistryLib
{
ChargeBalance createChargeBalance(BaseLib::ConfigTree const& config)
{
    auto const charge_balance =
        //! \ogs_file_param{prj__chemical_system__solution__charge_balance}
        config.getConfigParameterOptional<std::string>("charge_balance");

    if (!charge_balance)
    {
        return ChargeBalance::Unspecified;
    }
    if (*charge_balance == "pH")
    {
        return ChargeBalance::pH;
    }
    if (*charge_balance == "pe")
    {
        return ChargeBalance::pe;
    }

    OGS_FATAL(
        "Error in specifying the way in which charge balance is achieved.  "
        "Adjusting pH value or pe value are supported at this moment.");
}
}  // namespace ChemistryLib
