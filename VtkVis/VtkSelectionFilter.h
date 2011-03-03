/**
 * \file VtkSelectionFilter.h
 * 2010/02/09 KR Initial implementation
 *
 */


#ifndef VTKSELECTIONFILTER_H
#define VTKSELECTIONFILTER_H

// ** INCLUDES **
#include <vtkUnstructuredGridAlgorithm.h>
#include "VtkAlgorithmProperties.h"



class VtkSelectionFilter : public vtkUnstructuredGridAlgorithm, public VtkAlgorithmProperties
{

public:
	/// @brief Create new objects with New() because of VTKs object reference counting.
	static VtkSelectionFilter* New();

	vtkTypeRevisionMacro(VtkSelectionFilter, vtkUnstructuredGridAlgorithm);

	/// @brief Prints the mesh data to an output stream.
	void PrintSelf(ostream& os, vtkIndent indent);

	/// @brief Sets user properties.
 	void SetUserProperty(QString name, QVariant value)
 	{
		Q_UNUSED(name);
		Q_UNUSED(value);
 	}

	void SetSelectionArray(std::vector<double> selection, double threshold, bool ifSmaller=true);


protected:
	VtkSelectionFilter();
	~VtkSelectionFilter();

	/// @brief The filter logic.
	int RequestData(vtkInformation* request, 
		            vtkInformationVector** inputVector, 
					vtkInformationVector* outputVector);

private:
	std::vector<double> _selection;
	double _threshold;
	double _ifSmaller;
};

#endif // VTKSELECTIONFILTER_H
