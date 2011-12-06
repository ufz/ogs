/**
 * \file VtkCompositeGeoObjectFilter.h
 * 2011/12/02 KR Initial implementation
 */

#ifndef VTKCOMPOSITEGEOOBJECTFILTER_H
#define VTKCOMPOSITEGEOOBJECTFILTER_H

#include "VtkCompositeFilter.h"
#include "GeoType.h"

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

	void SetIndex(size_t idx, GEOLIB::GEOTYPE type) 
	{ 
		_index = idx; 
		_type = type;
	};

private:
	size_t _index;
	GEOLIB::GEOTYPE _type;
};

#endif // VTKCOMPOSITEGEOOBJECTFILTER_H
