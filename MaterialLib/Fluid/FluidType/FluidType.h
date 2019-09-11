/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

namespace FluidType
{
enum class Fluid_Type
{
    INCOMPRESSIBLE_FLUID,
    COMPRESSIBLE_FLUID,
    IDEAL_GAS
};

Fluid_Type strToFluidType(std::string const& s);

bool checkRequiredParams(Fluid_Type const& f_type,
                         double const& /*fluid_compressibility*/ beta_p,
                         double const& /*reference_temperature*/ T_ref,
                         double const& /*specific_gas_constant*/ R_s);

const char* getErrorMsg(Fluid_Type const& f_type);

}  // namespace FluidType
