/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
}  // namespace MeshLib

namespace ParameterLib
{
/// A parameter represented by a mesh property vector.
template <typename T>
struct MeshNodeParameter final : public Parameter<T>
{
    MeshNodeParameter(std::string const& name_,
                      MeshLib::Mesh const& mesh,
                      MeshLib::PropertyVector<T> const& property)
        : Parameter<T>(name_, &mesh), _property(property)
    {
    }

    bool isTimeDependent() const override { return false; }

    int getNumberOfGlobalComponents() const override
    {
        return _property.getNumberOfGlobalComponents();
    }

    std::vector<T> operator()(double const /*t*/,
                              SpatialPosition const& pos) const override
    {
        auto const n = pos.getNodeID();
        if (!n)
        {
            OGS_FATAL(
                "Trying to access a MeshNodeParameter but the node id is not "
                "specified.");
        }
        auto const num_comp = _property.getNumberOfGlobalComponents();
        std::vector<T> cache(num_comp);
        for (int c = 0; c < num_comp; ++c)
        {
            cache[c] = _property.getComponent(*n, c);
        }

        if (!this->_coordinate_system)
        {
            return cache;
        }

        return this->rotateWithCoordinateSystem(cache, pos);
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getNodalValuesOnElement(
        MeshLib::Element const& element, double const t) const override
    {
        auto const n_nodes = element.getNumberOfNodes();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(
            n_nodes, getNumberOfGlobalComponents());

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
};

std::unique_ptr<ParameterBase> createMeshNodeParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh);

}  // namespace ParameterLib
