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
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

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

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getNodalValuesOnElement(
        MeshLib::Element const& element, double const t) const override
    {
        auto const n_nodes = element.getNumberOfNodes();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(
            n_nodes, getNumberOfComponents());

        SpatialPosition x_position;
        auto const nodes = element.getNodes();
        for (unsigned i = 0; i < n_nodes; ++i)
        {
            x_position.setNodeID(nodes[i]->getID());
            auto const& values = this->operator()(t, x_position);
            result.row(i) =
                Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> const>(
                    values.data(), values.size());
        }

        return result;
    }

private:
    MeshLib::PropertyVector<T> const& _property;
    mutable std::vector<double> _cache;
};

std::unique_ptr<ParameterBase> createMeshNodeParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh);

}  // ProcessLib
