/**
 * \file VtkCompositeGeoObjectFilter.h
 * 2011/12/02 KR Initial implementation
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
	GEOLIB::GEOTYPE _type;
	vtkThreshold* _threshold;
};

#endif // VTKCOMPOSITEGEOOBJECTFILTER_H
