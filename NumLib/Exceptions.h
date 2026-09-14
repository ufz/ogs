// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <stdexcept>
#include <string>

namespace NumLib
{
/*! Thrown by a local assembler that cannot compute its contribution, e.g.,
 * because a local nonlinear solver did not converge.
 *
 * This is not a fatal error for the Picard and the Newton nonlinear solver of
 * NumLib::NonlinearSolver. Both catch it, abort the current nonlinear
 * iteration and report the time step as not converged, so that the time
 * stepping algorithm may repeat the step with a smaller time step size.
 *
 * \todo NumLib::PETScNonlinearSolver does not handle it. There the assembly
 * runs inside a callback called by PETSc's SNESSolve(), out of which a C++
 * exception cannot be thrown, so a process using the PETScSNES nonlinear
 * solver still ends the run.
 *
 * The exception must be thrown either on all MPI ranks or on none of them.
 * Ranks that continue while another one aborts run into the collective calls
 * of the linear solve, and the simulation deadlocks instead of failing.
 * Assemblers based on ProcessLib::AssemblyMixin get that behaviour from
 * BaseLib::MPI::allRanksThrowOrNone(); a hand-written assembler has to
 * establish it itself.
 */
struct AssemblyException : public std::runtime_error
{
    explicit AssemblyException(std::string const& reason)
        : std::runtime_error{"Error in process' assembly: " + reason} {};
};
}  // namespace NumLib
