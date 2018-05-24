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

    int getNumberOfComponents() const override
    {
        return _property.getNumberOfComponents();
    }

    std::vector<T> const& operator()(double const /*t*/,
                                     SpatialPosition const& pos) const override
    {
        auto const e = pos.getElementID();
        if (!e)
        {
            OGS_FATAL(
                "Trying to access a MeshElementParameter but the element id is "
                "not specified.");
        }
        auto const num_comp = _property.getNumberOfComponents();
        for (int c = 0; c < num_comp; ++c)
        {
            _cache[c] = _property.getComponent(*e, c);
        }
        return _cache;
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getNodalValuesOnElement(
        MeshLib::Element const& element, double const t) const override
    {
        auto const n_nodes = element.getNumberOfNodes();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(
            n_nodes, getNumberOfComponents());

        // Column vector of values, copied for each node.
        SpatialPosition x_position;
        x_position.setElementID(element.getID());
        auto const& values = this->operator()(t, x_position);
        auto const row_values =
            Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> const>(
                values.data(), values.size());
        for (unsigned i = 0; i < n_nodes; ++i)
        {
            result.row(i) = row_values;
        }
        return result;
    }

private:
    MeshLib::PropertyVector<T> const& _property;
    mutable std::vector<T> _cache;
};

std::unique_ptr<ParameterBase> createMeshElementParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh);

}  // ProcessLib
