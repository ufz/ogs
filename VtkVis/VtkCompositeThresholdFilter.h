/**
 * \file VtkCompositeThresholdFilter.h
 * 25/10/2010 LB Initial implementation
 */

#ifndef VTKCOMPOSITETHRESHOLDFILTER_H
#define VTKCOMPOSITETHRESHOLDFILTER_H

#include "VtkCompositeFilter.h"

class VtkCompositeThresholdFilter : public VtkCompositeFilter
{
public:
	VtkCompositeThresholdFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeThresholdFilter();
	
	virtual void init();

	virtual void SetUserProperty(QString name, QVariant value);

	void SetUserVectorProperty( QString name, QList<QVariant> values );

private:
	
};

#endif // VTKCOMPOSITETHRESHOLDFILTER_H
