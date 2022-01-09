/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <bitset>
#include <vector>

#include "BaseLib/Error.h"

#include "MeshLib/Node.h"

namespace FileIO
{
namespace Gocad
{
enum class FaceDirection : char
{
    U,
    V,
    W
};

class GocadNode : public MeshLib::Node
{
public:
    GocadNode(double const* const coords, std::size_t id,
              std::size_t layer_transition_idx)
        : Node(coords, id), _layer_transition_idx(layer_transition_idx)
    {
    }

    GocadNode() = delete;
    GocadNode(GocadNode const& src) = default;
    GocadNode(GocadNode&& src) = default;
    GocadNode& operator=(GocadNode&& rhs) = default;
    GocadNode& operator=(GocadNode const& rhs) = default;
    ~GocadNode() override = default;

    void setFaceSet(std::size_t face_set_number, std::size_t face_indicator)
    {
        _face_set_membership.set(face_set_number);
        switch (face_indicator)
        {
            case 0:
                _face_directions.emplace_back(face_set_number,
                                              FaceDirection::U);
                break;
            case 1:
                _face_directions.emplace_back(face_set_number,
                                              FaceDirection::V);
                break;
            case 2:
                _face_directions.emplace_back(face_set_number,
                                              FaceDirection::W);
                break;
            default:
                OGS_FATAL(
                    "GocadNode::setFaceSet(): unknown face indicator {:d}.",
                    face_indicator);
        }
    }

    /**
     * Checks if this GocadNode is in the face set with the number
     * face_set_number.
     * @param face_set_number the number of the face set
     * @return true/false
     */
    bool isMemberOfFaceSet(std::size_t face_set_number) const
    {
        return _face_set_membership[face_set_number];
    }

    bool isMemberOfAnyFaceSet() const { return _face_set_membership.any(); }

    void resetID(std::size_t id) { this->setID(id); }

    std::bitset<128> const& getFaceSetMembership() const
    {
        return _face_set_membership;
    }

    FaceDirection getFaceDirection(std::size_t const face_set_number) const
    {
        auto const it = std::find_if(
            _face_directions.begin(), _face_directions.end(),
            [&](auto const fi) { return fi.first == face_set_number; });
        if (it == _face_directions.end())
        {
            OGS_FATAL(
                "GocadNode {:d}: Could not found face indicator for face set "
                "{:d}",
                getID(), face_set_number);
        }
        return it->second;
    }

    std::size_t getLayerTransitionIndex() const
    {
        return _layer_transition_idx;
    }

protected:
    friend class GocadSplitNode;
    std::vector<std::pair<std::size_t, FaceDirection>> _face_directions;

private:
    std::bitset<128> _face_set_membership;
    std::size_t _layer_transition_idx;
};

bool operator<=(GocadNode const& n0, GocadNode const& n1);

class GocadSplitNode final : public GocadNode
{
public:
    GocadSplitNode(double const* const coords, std::size_t id,
                   std::array<std::size_t, 3> const& grid_coords,
                   std::bitset<8> const& affected_cells,
                   std::size_t layer_transition_idx)
        : GocadNode(coords, id, layer_transition_idx),
          _grid_coords(grid_coords),
          _affected_cells(affected_cells)
    {
    }

    std::array<std::size_t, 3> const& getGridCoords() const
    {
        return _grid_coords;
    }
    std::bitset<8> const& getAffectedCells() const
    {
        return _affected_cells;
    }
    void transmitFaceDirections(GocadNode const& gocad_node)
    {
        _face_directions = gocad_node._face_directions;
    }

private:
    std::array<std::size_t, 3> _grid_coords;
    std::bitset<8> const _affected_cells;
};

}  // end namespace Gocad
}  // end namespace FileIO

