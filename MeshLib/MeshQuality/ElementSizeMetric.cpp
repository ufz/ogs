/**
 * \file   ElementSizeMetric.cpp
 * \author Karsten Rink
 * \date   2011-03-17
 * \brief  Implementation of the ElementSizeMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    if (_mesh.getDimension() == 1)
        error_count = calc1dQuality();
    else if (_mesh.getDimension() == 2)
        error_count = calc2dQuality();
    else if (_mesh.getDimension() == 3)
        error_count = calc3dQuality();

    INFO ("ElementSizeMetric::calculateQuality() minimum: %f, max_volume: %f", _min, _max);
    if (error_count > 0)
        WARN ("Warning: %d elements with zero volume found.", error_count);
}

std::size_t ElementSizeMetric::calc1dQuality()
{
    const std::vector<MeshLib::Element*> &elements(_mesh.getElements());
    const std::size_t nElems(elements.size());
    std::size_t error_count(0);

    for (std::size_t k(0); k < nElems; k++)
    {
        double area(std::numeric_limits<double>::max());
        _element_quality_metric[k] = elements[k]->getContent();
        if (_element_quality_metric[k] < sqrt(fabs(std::numeric_limits<double>::epsilon())))
            error_count++;

        // update _min and _max values
        if (_min > area) _min = area;
        if (_max < area) _max = area;
    }
    return error_count;
}

std::size_t ElementSizeMetric::calc2dQuality()
{
    const std::vector<MeshLib::Element*> &elements(_mesh.getElements());
    const std::size_t nElems(elements.size());
    std::size_t error_count(0);

    for (std::size_t k(0); k < nElems; k++)
    {
        Element const& elem (*elements[k]);

        if (elem.getDimension() == 1)
        {
            _element_quality_metric[k] = 0.0;
            continue;
        }
        double const area = elem.getContent();
        if (area < sqrt(fabs(std::numeric_limits<double>::epsilon())))
            error_count++;

        // update _min and _max values
        if (_min > area) _min = area;
        if (_max < area) _max = area;
        _element_quality_metric[k] = area;
    }
    return error_count;
}

std::size_t ElementSizeMetric::calc3dQuality()
{
    const std::vector<MeshLib::Element*>& elements(_mesh.getElements());
    const std::size_t nElems(elements.size());
    std::size_t error_count(0);

    for (std::size_t k(0); k < nElems; k++)
    {
        Element const& elem (*elements[k]);
        if (elem.getDimension()<3)
        {
            _element_quality_metric[k] = 0.0;
            continue;
        }

        double const volume (elem.getContent());
        if (volume < sqrt(fabs(std::numeric_limits<double>::epsilon())))
            error_count++;

        if (_min > volume) _min = volume;
        if (_max < volume) _max = volume;
        _element_quality_metric[k] = volume;
    }
    return error_count;
}

} // end namespace MeshLib
