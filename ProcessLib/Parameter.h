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

#include <logog/include/logog.hpp>
#include <boost/optional.hpp>

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Elements/Element.h"
#include "SpatialPosition.h"

namespace MeshLib
{
template <typename T>
class PropertyVector;
}

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

/// Single, constant value parameter.
template <typename T>
struct ConstParameter final : public Parameter<T> {
    ConstParameter(T const& value) : _value{{value}} {}
    std::vector<T> const& getTuple(
        double const /*t*/, SpatialPosition const& /*pos*/) const override
    {
        return _value;
    }

private:
    std::vector<T> _value;
};

std::unique_ptr<ParameterBase> createConstParameter(BaseLib::ConfigTree const& config);

/// A parameter represented by a mesh property vector.
template <typename T>
struct MeshElementParameter final
    : public Parameter<T>
{
    MeshElementParameter(MeshLib::PropertyVector<T> const& property)
        : _property(property)
    {
    }

    std::vector<T> const& getTuple(
        double const /*t*/, SpatialPosition const& pos) const override
    {
        auto const e = pos.getElementID();
        assert(e);
        _cache.front() = _property[*e];
        return _cache;
    }

private:
    MeshLib::PropertyVector<T> const& _property;
    // TODO multi-component
    mutable std::vector<double> _cache = std::vector<double>(1);
};

std::unique_ptr<ParameterBase> createMeshPropertyParameter(BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh);

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PARAMETER_H_
