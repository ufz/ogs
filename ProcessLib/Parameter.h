/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PARAMETER_H_
#define PROCESS_LIB_PARAMETER_H_

#include <memory>
#include <vector>
#include "SpatialPosition.h"

namespace BaseLib
{
class ConfigTree;
}  // BaseLib

namespace MeshLib
{
class Mesh;
}  // MeshLib

namespace ProcessLib
{
/// Base class for all parameters, not an interface class. This avoids using of
/// void* when storing parameters and convenient destruction.
/// Its property name helps addressing the right parameter.
struct ParameterBase
{
    virtual ~ParameterBase() = default;

    std::string name;
};

template <typename T>
struct Parameter : public ParameterBase
{
    virtual ~Parameter() = default;

    // TODO number of components
    virtual std::vector<T> const& getTuple(
        double const t, SpatialPosition const& pos) const = 0;
};

std::unique_ptr<ParameterBase> createParameter(
    BaseLib::ConfigTree const& config,
    const std::vector<MeshLib::Mesh*>& meshes);

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PARAMETER_H_
