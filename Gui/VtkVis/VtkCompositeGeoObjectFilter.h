/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file VtkCompositeGeoObjectFilter.h
 *
 * Created on 2011-12-02 by Karsten Rink
 */

#ifndef VTKCOMPOSITEGEOOBJECTFILTER_H
#define VTKCOMPOSITEGEOOBJECTFILTER_H

#include "VtkCompositeFilter.h"
#include "GeoType.h"

class vtkThreshold;

/// @brief Hightlights a single GeoObject
class VtkCompositeGeoObjectFilter : public VtkCompositeFilter
{
public:
	VtkCompositeGeoObjectFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeGeoObjectFilter();

	virtual void init();

	/// @brief Sets user properties.
	void SetUserProperty(QString name, QVariant value)
	{
		Q_UNUSED(name);
		Q_UNUSED(value);
	}

	void SetIndex(size_t idx);

private:
	float GetInitialRadius() const;

	GeoLib::GEOTYPE _type;
	vtkThreshold* _threshold;
};

#endif // VTKCOMPOSITEGEOOBJECTFILTER_H
