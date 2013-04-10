/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-10
 * \brief  Definition of the VtkCompositeSelectionFilter class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKCOMPOSITESELECTIONFILTER_H
#define VTKCOMPOSITESELECTIONFILTER_H

#include "VtkCompositeFilter.h"

#include <vector>

class vtkThreshold;

class VtkColorLookupTable;

/// @brief This filter colors the input by the points z-value.
class VtkCompositeSelectionFilter : public VtkCompositeFilter
{
public:
	VtkCompositeSelectionFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeSelectionFilter() {}

	virtual void init();

	void setSelectionArray(const std::string &selection_name, bool is_element_array = true, const std::vector<double> &selection = std::vector<double>());

	virtual void SetUserVectorProperty(QString name, QList<QVariant> values);

private:
	/// Returns a colour lookup table optimised for quality measures
	VtkColorLookupTable* GetLookupTable();

	vtkThreshold* _threshold;
	std::string _selection_name;
	std::vector<double> _selection;
	bool _is_element_array;
};

#endif // VTKCOMPOSITESELECTIONFILTER_H
