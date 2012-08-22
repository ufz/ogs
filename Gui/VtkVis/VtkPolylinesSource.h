/**
 * \file VtkPolylinesSource.h
 * 2/2/2010 LB Initial implementation
 *
 */

#ifndef VTKPOLYLINESSOURCE_H
#define VTKPOLYLINESSOURCE_H

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include <vtkPolyDataAlgorithm.h>

#include "GEOObjects.h"

/**
 * \brief VtkPolylinesSource is a VTK source object for the visualisation of
 * polyline data. As a vtkPolyDataAlgorithm it outputs polygonal data.
 */
class VtkPolylinesSource : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
	/// Create new objects with New() because of VTKs object reference counting.
	static VtkPolylinesSource* New();

	vtkTypeRevisionMacro(VtkPolylinesSource,vtkPolyDataAlgorithm);

	/// Sets the polyline vector.
	void setPolylines(const std::vector<GeoLib::Polyline*>* polylines) { _polylines = polylines; }

	/// Prints its data on a stream.
	void PrintSelf(ostream& os, vtkIndent indent);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	VtkPolylinesSource();
	~VtkPolylinesSource();

	/// Computes the polygonal data object.
	int RequestData(vtkInformation* request,
	                vtkInformationVector** inputVector,
	                vtkInformationVector* outputVector);

	int RequestInformation(vtkInformation* request,
	                       vtkInformationVector** inputVector,
	                       vtkInformationVector* outputVector);

	/// The polylines to visualize.
	const std::vector<GeoLib::Polyline*>* _polylines;

private:
};

#endif // VTKPOLYLINESSOURCE_H
