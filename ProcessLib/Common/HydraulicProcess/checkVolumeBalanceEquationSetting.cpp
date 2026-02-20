// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "checkVolumeBalanceEquationSetting.h"

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/Constant.h"
#include "MaterialLib/MPL/PropertyType.h"

namespace ProcessLib::Common::HydraulicProcess
{

void checkVolumeBalanceEquationSetting(
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    // Check whether the fluid phase density is constant for all media.
    for (auto const& medium : media_map.media())
    {
        // auto const& medium = *media_map.getMedium(element_id);
        auto const& fluid_phase_density = MaterialPropertyLib::fluidPhase(
            *medium)[MaterialPropertyLib::PropertyType::density];
        if (typeid(fluid_phase_density) !=
            typeid(MaterialPropertyLib::Constant))
        {
            OGS_FATAL(
                "Since `equation_balance_type` is set to `volume`,the "
                "phase density type must be `Constant`. Note: by "
                "default, `equation_balance_type` is set to `volume`.");
        }
    }
}
}  // namespace ProcessLib::Common::HydraulicProcess
