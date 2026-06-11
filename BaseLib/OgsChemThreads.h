// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace BaseLib
{
/// Returns the number of threads to use for parallel chemistry calculations.
/// Priority: OGS_CHEM_THREADS env var, then OGS_ASM_THREADS env var, then 1.
int getNumberOfChemistryThreads();
}  // namespace BaseLib
