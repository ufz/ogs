/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_INITIAL_CONDITION_H_
#define PROCESS_LIB_INITIAL_CONDITION_H_

#include <cassert>
#include "BaseLib/ConfigTree.h"
#include "MeshLib/PropertyVector.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
template <typename>
class PropertyVector;
class Mesh;
}

namespace ProcessLib
{
/// The InitialCondition is a base class for spatial distributions of values
/// defined on mesh nodes.
class InitialCondition
{
public:
    virtual ~InitialCondition() = default;
    virtual double getValue(std::size_t const node_id,
                            int const component_id) const = 0;
};

/// Uniform value initial condition
class UniformInitialCondition : public InitialCondition
{
public:
    explicit UniformInitialCondition(std::vector<double> const& values)
        : _values(values)
    {
    }
    /// Returns a value for given node and component.
    virtual double getValue(std::size_t const /*node_id*/,
                            int const component_id) const override
    {
        return _values[component_id];
    }

private:
    std::vector<double> const _values;
};

/// Construct a UniformInitialCondition from configuration.
/// The initial condition will expect a correct number of components in the
/// configuration, which should be the same as in the corresponding process
/// variable.
std::unique_ptr<InitialCondition> createUniformInitialCondition(
    BaseLib::ConfigTree const& config, int const n_components);

/// Distribution of values given by a mesh property defined on nodes.
class MeshPropertyInitialCondition : public InitialCondition
{
public:
    MeshPropertyInitialCondition(
        MeshLib::PropertyVector<double> const& property)
        : _property(property)
    {
        assert(_property.getMeshItemType() == MeshLib::MeshItemType::Node);
    }

    virtual double getValue(std::size_t const node_id,
                            int const component_id) const override
    {
        return _property.getComponent(node_id, component_id);
    }

private:
    MeshLib::PropertyVector<double> const& _property;
};

/// Construct a MeshPropertyInitialCondition from configuration.
/// The initial condition will expect a correct number of components in the
/// configuration, which should be the same as in the corresponding process
/// variable.
std::unique_ptr<InitialCondition> createMeshPropertyInitialCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh,
    int const n_components);

}  // namespace ProcessLib

#endif  // PROCESS_LIB_INITIAL_CONDITION_H_
