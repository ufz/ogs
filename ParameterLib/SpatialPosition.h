/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/optional.hpp>
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
        boost::optional<std::size_t> const& node_id,
        boost::optional<std::size_t> const& element_id,
        boost::optional<unsigned> const& integration_point,
        boost::optional<MathLib::TemplatePoint<double, 3>> const& coordinates)
        : node_id_(node_id),
          element_id_(element_id),
          integration_point_(integration_point),
          coordinates_(coordinates)
    {
    }

    boost::optional<std::size_t> getNodeID() const { return node_id_; }
    boost::optional<std::size_t> getElementID() const { return element_id_; }
    boost::optional<unsigned> getIntegrationPoint() const
    {
        return integration_point_;
    }
    boost::optional<MathLib::TemplatePoint<double, 3>> const& getCoordinates()
        const
    {
        return coordinates_;
    }

    void setNodeID(std::size_t node_id)
    {
        clear();
        node_id_ = node_id;
    }

    void setElementID(std::size_t element_id)
    {
        clear();
        element_id_ = element_id;
    }

    void setIntegrationPoint(unsigned integration_point)
    {
        assert(element_id_);
        integration_point_ = integration_point;
    }

    void setCoordinates(MathLib::TemplatePoint<double, 3> const& coordinates)
    {
        coordinates_ = coordinates;
    }

    void setAll(
        boost::optional<std::size_t> const& node_id,
        boost::optional<std::size_t> const& element_id,
        boost::optional<unsigned> const& integration_point,
        boost::optional<MathLib::TemplatePoint<double, 3>> const& coordinates)
    {
        node_id_ = node_id;
        element_id_ = element_id;
        integration_point_ = integration_point;
        coordinates_ = coordinates;
    }

    void clear()
    {
        node_id_ = boost::none;
        element_id_ = boost::none;
        integration_point_ = boost::none;
        coordinates_ = boost::none;
    }

private:
    boost::optional<std::size_t> node_id_;
    boost::optional<std::size_t> element_id_;
    boost::optional<unsigned> integration_point_;
    boost::optional<MathLib::TemplatePoint<double, 3>> coordinates_;
};

}  // namespace ParameterLib
