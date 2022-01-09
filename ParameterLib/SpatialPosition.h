/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <optional>

#include "MathLib/TemplatePoint.h"

namespace ParameterLib
{
//! Represents a position in space which can be either one of
//! a node, an element, an integration point or a cartesian coordinates triple.
//!
//! The setters of this class make sure that only compatible information can be
//! stored at the same time; e.g., it is not possible to specify an element ID
//! and a node ID at the same time (the setAll() method being an exception to
//! that rule).
class SpatialPosition
{
public:
    SpatialPosition() = default;

    SpatialPosition(
        std::optional<std::size_t> const& node_id,
        std::optional<std::size_t> const& element_id,
        std::optional<unsigned> const& integration_point,
        std::optional<MathLib::TemplatePoint<double, 3>> const& coordinates)
        : _node_id(node_id),
          _element_id(element_id),
          _integration_point(integration_point),
          _coordinates(coordinates)
    {
    }

    std::optional<std::size_t> getNodeID() const { return _node_id; }
    std::optional<std::size_t> getElementID() const { return _element_id; }
    std::optional<unsigned> getIntegrationPoint() const
    {
        return _integration_point;
    }
    std::optional<MathLib::TemplatePoint<double, 3>> const& getCoordinates()
        const
    {
        return _coordinates;
    }

    void setNodeID(std::size_t node_id)
    {
        clear();
        _node_id = node_id;
    }

    void setElementID(std::size_t element_id)
    {
        clear();
        _element_id = element_id;
    }

    void setIntegrationPoint(unsigned integration_point)
    {
        assert(_element_id);
        _integration_point = integration_point;
    }

    void setCoordinates(MathLib::TemplatePoint<double, 3> const& coordinates)
    {
        _coordinates = coordinates;
    }

    void setAll(
        std::optional<std::size_t> const& node_id,
        std::optional<std::size_t> const& element_id,
        std::optional<unsigned> const& integration_point,
        std::optional<MathLib::TemplatePoint<double, 3>> const& coordinates)
    {
        _node_id = node_id;
        _element_id = element_id;
        _integration_point = integration_point;
        _coordinates = coordinates;
    }

    void clear()
    {
        _node_id = std::nullopt;
        _element_id = std::nullopt;
        _integration_point = std::nullopt;
        _coordinates = std::nullopt;
    }

private:
    std::optional<std::size_t> _node_id;
    std::optional<std::size_t> _element_id;
    std::optional<unsigned> _integration_point;
    std::optional<MathLib::TemplatePoint<double, 3>> _coordinates;
};

}  // namespace ParameterLib
