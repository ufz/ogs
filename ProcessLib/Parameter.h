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

/// A parameter is representing a value or function of any type.
/// The ReturnType can represent an underlying type of an aggregate type like
/// tuple or matrix (\see tupleSize()).
/// The total number of stored tuples is provided.
template <typename ReturnType, typename... Args>
struct Parameter : public ParameterBase
{
    virtual ~Parameter() = default;

    virtual ReturnType operator()(Args&&... args) const = 0;
};

/// Single, constant value parameter.
template <typename ReturnType>
struct ConstParameter final
    : public Parameter<ReturnType, MeshLib::Element const&>
{
    ConstParameter(ReturnType value) : _value(value)
    {
    }

    ReturnType operator()(MeshLib::Element const&) const override
    {
        return _value;
    }

private:
    ReturnType _value;
};

std::unique_ptr<ParameterBase> createConstParameter(BaseLib::ConfigTree const& config);

/// A parameter represented by a mesh property vector.
template <typename ReturnType>
struct MeshPropertyParameter final
    : public Parameter<ReturnType, MeshLib::Element const&>
{
    MeshPropertyParameter(MeshLib::PropertyVector<ReturnType> const& property)
        : _property(property)
    {
    }

    ReturnType operator()(MeshLib::Element const& e) const override
    {
        return _property[e.getID()];
    }

private:
    MeshLib::PropertyVector<ReturnType> const& _property;
};

std::unique_ptr<ParameterBase> createMeshPropertyParameter(BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh);

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PARAMETER_H_
