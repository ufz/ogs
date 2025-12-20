// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
