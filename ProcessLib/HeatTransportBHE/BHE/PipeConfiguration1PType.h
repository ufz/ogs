/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Pipe.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct PipeConfiguration1PType
{
    Pipe const single_pipe;

    double const longitudinal_dispersion_length;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
