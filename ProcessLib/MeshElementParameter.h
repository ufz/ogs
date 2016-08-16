/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_MESHELEMENTPARAMETER_H
#define PROCESSLIB_MESHELEMENTPARAMETER_H

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
    MeshElementParameter(MeshLib::PropertyVector<T> const& property)
        : _property(property)
    {
    }

    std::vector<T> const& getTuple(double const /*t*/,
                                   SpatialPosition const& pos) const override
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

std::unique_ptr<ParameterBase> createMeshElementParameter(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh);

}  // ProcessLib

#endif  // PROCESSLIB_MESHELEMENTPARAMETER_H
