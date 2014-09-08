/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-10
 * \brief  Definition of the VtkCompositeSelectionFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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

/// @brief This filter selects/thresholds elements based on the selected array.
class VtkCompositeElementSelectionFilter : public VtkCompositeFilter
{
public:
	VtkCompositeElementSelectionFilter(vtkAlgorithm* inputAlgorithm);
	virtual ~VtkCompositeElementSelectionFilter() {}

	virtual void init();

	void setSelectionArray(const std::string &selection_name, const std::vector<double> &selection = std::vector<double>());

	virtual void SetUserVectorProperty(QString name, QList<QVariant> values);

private:
	/// Returns a colour lookup table optimised for quality measures
	VtkColorLookupTable* GetLookupTable();

	std::string _selection_name;
	std::vector<double> _selection;
};

#endif // VTKCOMPOSITESELECTIONFILTER_H
