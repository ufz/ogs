/**
 * \file   AreaMetric.cpp
 * \author Karsten Rink
 * \date   2011-03-17
 * \brief  Implementation of the AreaMetric class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include "AreaMetric.h"
#include "MathTools.h"

namespace MeshLib
{
AreaMetric::AreaMetric(Mesh const* const mesh)
	: ElementQualityMetric(mesh)
{}

void AreaMetric::calculateQuality()
{
	const std::vector<MeshLib::Element*> &elements(_mesh->getElements());

	const size_t nElems(elements.size());
	for (size_t k(0); k < nElems; k++) 
	{
		double area(std::numeric_limits<double>::max());
		const Element* elem (elements[k]);

		if (elem->getDimension() == 1)
		{
			_element_quality_metric[k] = -1.0;
			continue;
		}
		else if (elem->getDimension() == 2)
		{		
			area = elem->getContent();
			if (area < sqrt(fabs(std::numeric_limits<double>::epsilon()))) errorMsg(elem, k);
		} 
		else {
			size_t nFaces(elem->getNFaces());

			for (size_t i = 0; i < nFaces; i++) 
			{
				const double sub_area (elem->getFace(i)->getContent());

				if (sub_area < sqrt(fabs(std::numeric_limits<double>::epsilon())))
					errorMsg(elem, k);
				if (sub_area < area) area = sub_area;
			}
		}
		// update _min and _max values
		if (_min > area) _min = area;
		if (_max < area) _max = area;
		_element_quality_metric[k] = area;
	}
}

} // end namespace MeshLib
