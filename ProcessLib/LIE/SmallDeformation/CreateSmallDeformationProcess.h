/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_SMALLDEFORMATION_CREATESMALLDEFORMATIONPROCESS_H_
#define PROCESSLIB_LIE_SMALLDEFORMATION_CREATESMALLDEFORMATIONPROCESS_H_

#include <memory>
#include <vector>

#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{

template <int DisplacementDim>
std::unique_ptr<Process> createSmallDeformationProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib

#endif  // PROCESSLIB_LIE_SMALLDEFORMATION_CREATESMALLDEFORMATIONPROCESS_H_
