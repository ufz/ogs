/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_CREATESMALLDEFORMATIONPROCESS_H_
#define PROCESS_LIB_CREATESMALLDEFORMATIONPROCESS_H_

#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace SmallDeformation
{

template <int DisplacementDim>
std::unique_ptr<Process>
createSmallDeformationProcess(
    MeshLib::Mesh& mesh,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}  // namespace SmallDeformation
}  // namespace ProcessLib

#endif  // PROCESS_LIB_CREATESMALLDEFORMATIONPROCESS_H_
