/**
 * \file
 * \author Karsten Rink
 * \date   2011-03-17
 * \brief  Implementation of the ElementSizeMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementSizeMetric.h"

#include <limits>

namespace MeshLib
{
ElementSizeMetric::ElementSizeMetric(Mesh const& mesh)
: ElementQualityMetric(mesh)
{}

void ElementSizeMetric::calculateQuality()
{
    std::size_t error_count(0);
    if (mesh_.getDimension() == 1)
    {
        error_count = calc1dQuality();
    }
    else if (mesh_.getDimension() == 2)
    {
        error_count = calc2dQuality();
    }
    else if (mesh_.getDimension() == 3)
    {
        error_count = calc3dQuality();
    }

    INFO(
        "ElementSizeMetric::calculateQuality() minimum: {:f}, max_volume: {:f}",
        min_,
        max_);
    if (error_count > 0)
    {
        WARN("Warning: {:d} elements with zero volume found.", error_count);
    }
}

std::size_t ElementSizeMetric::calc1dQuality()
{
    const std::vector<MeshLib::Element*> &elements(mesh_.getElements());
    const std::size_t nElems(elements.size());
    std::size_t error_count(0);

    for (std::size_t k(0); k < nElems; k++)
    {
        double area(std::numeric_limits<double>::max());
        element_quality_metric_[k] = elements[k]->getContent();
        if (element_quality_metric_[k] <
            sqrt(fabs(std::numeric_limits<double>::epsilon())))
        {
            error_count++;
        }

        // update min_ and max_ values
        if (min_ > area)
        {
            min_ = area;
        }
        if (max_ < area)
        {
            max_ = area;
        }
    }
    return error_count;
}

std::size_t ElementSizeMetric::calc2dQuality()
{
    const std::vector<MeshLib::Element*> &elements(mesh_.getElements());
    const std::size_t nElems(elements.size());
    std::size_t error_count(0);

    for (std::size_t k(0); k < nElems; k++)
    {
        Element const& elem (*elements[k]);

        if (elem.getDimension() == 1)
        {
            element_quality_metric_[k] = 0.0;
            continue;
        }
        double const area = elem.getContent();
        if (area < sqrt(fabs(std::numeric_limits<double>::epsilon())))
        {
            error_count++;
        }

        // update min_ and max_ values
        if (min_ > area)
        {
            min_ = area;
        }
        if (max_ < area)
        {
            max_ = area;
        }
        element_quality_metric_[k] = area;
    }
    return error_count;
}

std::size_t ElementSizeMetric::calc3dQuality()
{
    const std::vector<MeshLib::Element*>& elements(mesh_.getElements());
    const std::size_t nElems(elements.size());
    std::size_t error_count(0);

    for (std::size_t k(0); k < nElems; k++)
    {
        Element const& elem (*elements[k]);
        if (elem.getDimension()<3)
        {
            element_quality_metric_[k] = 0.0;
            continue;
        }

        double const volume (elem.getContent());
        if (volume < sqrt(fabs(std::numeric_limits<double>::epsilon())))
        {
            error_count++;
        }

        if (min_ > volume)
        {
            min_ = volume;
        }
        if (max_ < volume)
        {
            max_ = volume;
        }
        element_quality_metric_[k] = volume;
    }
    return error_count;
}

} // end namespace MeshLib
