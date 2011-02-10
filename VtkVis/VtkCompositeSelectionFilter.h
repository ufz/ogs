/**
 * \file VtkCompositeSelectionFilter.h
 * 2011/02/10 KR Initial implementation
 */

#ifndef VTKCOMPOSITESELECTIONFILTER_H
#define VTKCOMPOSITESELECTIONFILTER_H

#include "VtkCompositeFilter.h"

/// @brief This filter colors the input by the points z-value.
class VtkCompositeSelectionFilter : public VtkCompositeFilter
{
public:
	VtkCompositeSelectionFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeSelectionFilter() {};

	virtual void init();

	void setSelectionArray(std::vector<double> selection) { _selection = selection; init(); };

	virtual void SetUserProperty(QString name, QVariant value);

private:
	std::vector<double> _selection;
	
};

#endif // VTKCOMPOSITESELECTIONFILTER_H
