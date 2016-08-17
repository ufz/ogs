/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SPATIALPOSITION_H
#define PROCESSLIB_SPATIALPOSITION_H

#include <boost/optional.hpp>
#include "MathLib/TemplatePoint.h"

namespace ProcessLib
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
    boost::optional<std::size_t> getNodeID() const { return _node_id; }
    boost::optional<std::size_t> getElementID() const { return _element_id; }
    boost::optional<unsigned> getIntegrationPoint() const
    {
        return _integration_point;
    }
    boost::optional<MathLib::TemplatePoint<double, 3>> const& getCoordinates()
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
        clear();
        _coordinates = coordinates;
    }

    void setAll(
        boost::optional<std::size_t> const& node_id,
        boost::optional<std::size_t> const& element_id,
        boost::optional<unsigned> const& integration_point,
        boost::optional<MathLib::TemplatePoint<double, 3>> const& coordinates)
    {
        _node_id = node_id;
        _element_id = element_id;
        _integration_point = integration_point;
        _coordinates = coordinates;
    }

    void clear()
    {
        _node_id = boost::none;
        _element_id = boost::none;
        _integration_point = boost::none;
        _coordinates = boost::none;
    }

private:
    boost::optional<std::size_t> _node_id;
    boost::optional<std::size_t> _element_id;
    boost::optional<unsigned> _integration_point;
    boost::optional<MathLib::TemplatePoint<double, 3>> _coordinates;
};

}  // ProcessLib

#endif  // PROCESSLIB_SPATIALPOSITION_H
