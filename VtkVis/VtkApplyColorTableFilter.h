/**
 * \file VtkApplyColorTableFilter.h
 * 21/10/2010 KR Initial implementation
 *
 */


#ifndef VTKAPPLYCOLORTABLEFILTER_H
#define VTKAPPLYCOLORTABLEFILTER_H

// ** INCLUDES **
#include <vtkPolyDataAlgorithm.h>
#include "VtkAlgorithmProperties.h"
#include <vtkLookupTable.h>

/**
 * \brief Applying a color table to a vtk object.
 */
class VtkApplyColorTableFilter : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{

public:
	/// Create new objects with New() because of VTKs object reference counting.
	static VtkApplyColorTableFilter* New();

	vtkTypeRevisionMacro(VtkApplyColorTableFilter, vtkPolyDataAlgorithm);

	/// @brief Prints information about itself.
	void PrintSelf(ostream& os, vtkIndent indent);

	/// Returns the underlying colour look up table object.
	vtkGetObjectMacro(ColorLookupTable,vtkLookupTable);

	/// Sets the underlying colour look up table object.
	virtual void SetColorLookupTable(vtkLookupTable*);

	/// @brief Colors the whole cell instead of interpolating colors between
	/// the points in a cell.
	ogsUserPropertyMacro(ColorsOnCells,bool);

	//void SetColorsOnPointsOn() { _colorsOnPoints = true; };
	//void SetColorsOnCellsOn() { _colorsOnPoints = false; };

	virtual unsigned long GetMTime();

	virtual void SetUserProperty(QString name, QVariant value)
	{
		if (name.compare("ColorsOnCells") == 0)
			SetColorsOnCells(value.toBool());
	}

protected:
	VtkApplyColorTableFilter();
	~VtkApplyColorTableFilter();

	/// Computes the unstructured grid data object.
	int RequestData(vtkInformation* request, 
		            vtkInformationVector** inputVector, 
					vtkInformationVector* outputVector);


private:
	vtkLookupTable* ColorLookupTable;
	//bool _colorsOnPoints;
};

#endif // VTKAPPLYCOLORTABLEFILTER_H
