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

#ifndef AREAMETRIC_H_
#define AREAMETRIC_H_

#include "ElementQualityMetric.h"

namespace MeshLib
{

/** 
 * Calculates the quality of mesh elements based on the area of element or its faces
 */
class AreaMetric : public ElementQualityMetric
{
public:
	AreaMetric(Mesh const* const mesh);
	virtual ~AreaMetric() {}

	virtual void calculateQuality ();
};
}

#endif /* AREAMETRIC_H_ */
