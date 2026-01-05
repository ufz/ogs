// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <bitset>
#include <optional>

#include "MathLib/Point3d.h"

namespace ParameterLib
{
//! Represents a position in space which can be either one of
//! a node, an element, an integration point or a cartesian coordinates triple.
//!
//! The setters of this class make sure that only compatible information can be
//! stored at the same time; e.g., it is not possible to specify an element ID
//! and a node ID at the same time (the constructor being an exception to
//! that rule).
class SpatialPosition
{
public:
    SpatialPosition() = default;

    SpatialPosition(std::optional<std::size_t> const& node_id,
                    std::optional<std::size_t> const& element_id,
                    std::optional<MathLib::Point3d> const& coordinates)
    {
        if (node_id)
        {
            _node_id = *node_id;
            flags.set(node_bit);
        }

        if (element_id)
        {
            _element_id = *element_id;
            flags.set(element_bit);
        }
        if (coordinates)
        {
            _coordinates = *coordinates;
            flags.set(coordinates_bit);
        }
    }

    std::optional<std::size_t> getNodeID() const
    {
        return flags[node_bit] ? std::make_optional(_node_id) : std::nullopt;
    }
    std::optional<std::size_t> getElementID() const
    {
        return flags[element_bit] ? std::make_optional(_element_id)
                                  : std::nullopt;
    }
    std::optional<MathLib::Point3d> const getCoordinates() const
    {
        return flags[coordinates_bit] ? std::make_optional(_coordinates)
                                      : std::nullopt;
    }

    void setNodeID(std::size_t node_id)
    {
        flags.reset();
        flags.set(node_bit);
        _node_id = node_id;
    }

    void setElementID(std::size_t element_id)
    {
        flags.reset();
        flags.set(element_bit);
        _element_id = element_id;
    }

    void setCoordinates(MathLib::Point3d const& coordinates)
    {
        _coordinates = coordinates;
        flags.set(coordinates_bit);
    }

private:
    std::size_t _node_id = 0;
    std::size_t _element_id = 0;
    MathLib::Point3d _coordinates{};

    std::bitset<3> flags{};
    static constexpr std::size_t node_bit = 0;
    static constexpr std::size_t element_bit = 1;
    static constexpr std::size_t coordinates_bit = 2;
};

}  // namespace ParameterLib
