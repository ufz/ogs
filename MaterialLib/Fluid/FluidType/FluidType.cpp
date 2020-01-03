/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FluidType.h"
#include "BaseLib/Error.h"
#include <math.h>

namespace FluidType
{

Fluid_Type strToFluidType(std::string const& s)
{
    if (s == "incompressible_fluid")
    {
        return Fluid_Type::INCOMPRESSIBLE_FLUID;
    }
    if (s == "compressible_fluid")
    {
        return Fluid_Type::COMPRESSIBLE_FLUID;
    }
    if (s == "ideal_gas")
    {
        return Fluid_Type::IDEAL_GAS;
    }
    OGS_FATAL("This fluid type is unavailable. The available types are \n"
               "incompressible_fluid, compressible_fluid and ideal_gas.\n");
}

bool checkRequiredParams(Fluid_Type const& f_type,
                         double const& /*fluid_compressibility*/ beta_p,
                         double const& /*reference_temperature*/ T_ref,
                         double const& /*specific_gas_constant*/ R_s)
{
    return !((f_type == Fluid_Type::COMPRESSIBLE_FLUID && isnan(beta_p)) ||
             (f_type == Fluid_Type::IDEAL_GAS &&
              (isnan(T_ref) || isnan(R_s))));
}

const char* getErrorMsg(Fluid_Type const& f_type)
{
    if (f_type == Fluid_Type::COMPRESSIBLE_FLUID)
    {
        return "The fluid type compressible_fluid requires the parameter\n"
               "fluid_compressibility";
    }
    if (f_type == Fluid_Type::IDEAL_GAS)
    {
        return "The fluid type ideal_gas requires the parameters\n"
               "reference_temperature and specific_gas_constant.";
    }
    return "The required parameters for this fluid type are missing.";
}

}  // namespace FluidType
