/**
 * \file   VolumeMetric.h
 * \author Thomas Fischer
 * \date   2011-03-03
 * \brief  Definition of the VolumeMetric class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VOLUMEMETRIC_H_
#define VOLUMEMETRIC_H_

#include "ElementQualityMetric.h"

namespace MeshLib
{

/** 
 * Calculates the quality of mesh elements based on the volume of elements
 */
class VolumeMetric : public ElementQualityMetric
{
public:
	VolumeMetric(Mesh const& mesh);
	virtual ~VolumeMetric() {}

	virtual void calculateQuality ();
};
}

#endif /* VOLUMEMETRIC_H_ */
