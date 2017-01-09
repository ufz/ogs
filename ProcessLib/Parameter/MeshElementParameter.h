/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Parameter.h"

namespace MeshLib
{
template <typename T>
class PropertyVector;
}  // MeshLib

namespace ProcessLib
{
/// A parameter represented by a mesh property vector.
template <typename T>
struct MeshElementParameter final : public Parameter<T> {
    MeshElementParameter(std::string const& name_,
                         MeshLib::PropertyVector<T> const& property)
        : Parameter<T>(name_),
          _property(property),
          _cache(_property.getNumberOfComponents())
    {
    }

    bool isTimeDependent() const override { return false; }

    unsigned getNumberOfComponents() const override
    {
        return _property.getNumberOfComponents();
    }

    std::vector<T> const& operator()(double const /*t*/,
                                     SpatialPosition const& pos) const override
    {
        auto const e = pos.getElementID();
        assert(e);
        auto const num_comp = _property.getNumberOfComponents();
        for (std::size_t c=0; c<num_comp; ++c) {
            _cache[c] = _property.getComponent(*e, c);
        }
        return _cache;
    }

private:
    MeshLib::PropertyVector<T> const& _property;
    mutable std::vector<double> _cache;
};

std::unique_ptr<ParameterBase> createMeshElementParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh);

}  // ProcessLib
