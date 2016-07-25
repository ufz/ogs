/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_CREATE_HEATTRANSPORTPROCESS_H_
#define PROCESS_LIB_CREATE_HEATTRANSPORTPROCESS_H_

#include <memory>
#include "ProcessLib/Process.h"


namespace ProcessLib
{
namespace HeatTransport
{

std::unique_ptr<Process>
createHeatTransportProcess(
    MeshLib::Mesh& mesh,
    Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Process::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}   // namespace HeatTransport
}   // namespace ProcessLib

#endif  // PROCESS_LIB_CREATE_HEATTRANSPORTPROCESS_H_
