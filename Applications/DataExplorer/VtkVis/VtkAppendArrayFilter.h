/**
 * \file
 * \author Karsten Rink
 * \date   2010-02-09
 * \brief  Definition of the VtkAppendArrayFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKAPPENDARRAYFILTER_H
#define VTKAPPENDARRAYFILTER_H

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"

#include <vtkUnstructuredGridAlgorithm.h>

#include <vector>

/**
 * \brief
 */
class VtkAppendArrayFilter : public vtkUnstructuredGridAlgorithm, public VtkAlgorithmProperties
{
public:
	/// @brief Create new objects with New() because of VTKs object reference counting.
	static VtkAppendArrayFilter* New();

	vtkTypeMacro(VtkAppendArrayFilter, vtkUnstructuredGridAlgorithm);

	/// @brief Prints the mesh data to an output stream.
	void PrintSelf(ostream& os, vtkIndent indent);

	/// @brief Sets user properties.
	void SetUserProperty(QString name, QVariant value)
	{
		Q_UNUSED(name);
		Q_UNUSED(value);
	}

	void SetArray(const std::string &array_name, const std::vector<double> &selection);

protected:
	VtkAppendArrayFilter();
	~VtkAppendArrayFilter();

	/// @brief The filter logic.
	int RequestData(vtkInformation* request,
	                vtkInformationVector** inputVector,
	                vtkInformationVector* outputVector);

private:
	std::vector<double> _array;
	std::string _array_name;
};

#endif // VTKAPPENDARRAYFILTER_H
