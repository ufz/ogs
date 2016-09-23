/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_MATERIALID_PARAMETER_H_
#define PROCESSLIB_MATERIALID_PARAMETER_H_

#include "Parameter.h"

#include "MeshLib/PropertyVector.h"

namespace MeshLib
{
template <typename T>
class PropertyVector;
}  // MeshLib

namespace ProcessLib
{

std::unique_ptr<ParameterBase> createMaterialIDParameter(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh);

}  // ProcessLib

#endif // PROCESSLIB_MATERIALID_PARAMETER_H_
