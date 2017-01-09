/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_CREATE_TESPROCESS_H_
#define PROCESS_LIB_CREATE_TESPROCESS_H_

#include <memory>
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace TES
{
std::unique_ptr<Process> createTESProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& /*parameters*/,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace TES
}  // namespace ProcessLib

#endif  // PROCESS_LIB_CREATE_TESPROCESS_H_
