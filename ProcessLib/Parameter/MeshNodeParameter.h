/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Parameter.h"

#include "BaseLib/Error.h"

namespace MeshLib
{
template <typename T>
class PropertyVector;
}  // MeshLib

namespace ProcessLib
{
/// A parameter represented by a mesh property vector.
template <typename T>
struct MeshNodeParameter final : public Parameter<T> {
    MeshNodeParameter(std::string const& name_,
                      MeshLib::PropertyVector<T> const& property)
        : Parameter<T>(name_),
          _property(property),
          _cache(_property.getNumberOfComponents())
    {
    }

    bool isTimeDependent() const override { return false; }

    int getNumberOfComponents() const override
    {
        return _property.getNumberOfComponents();
    }

    std::vector<T> const& operator()(double const /*t*/,
                                     SpatialPosition const& pos) const override
    {
        auto const n = pos.getNodeID();
        if (!n)
        {
            OGS_FATAL(
                "Trying to access a MeshNodeParameter but the node id is not "
                "specified.");
        }
        auto const num_comp = _property.getNumberOfComponents();
        for (int c = 0; c < num_comp; ++c)
        {
            _cache[c] = _property.getComponent(*n, c);
        }
        return _cache;
    }

private:
    MeshLib::PropertyVector<T> const& _property;
    mutable std::vector<double> _cache;
};

std::unique_ptr<ParameterBase> createMeshNodeParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh);

}  // ProcessLib
