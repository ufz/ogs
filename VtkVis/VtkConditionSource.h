/**
 * \file VtkConditionSource.h
 * 2011/03/02 KR Initial implementation
 *
 */


#ifndef VTKCONDITIONSOURCE_H
#define VTKCONDITIONSOURCE_H

// ** INCLUDES **
#include <vtkPolyDataAlgorithm.h>
#include "VtkAlgorithmProperties.h"

#include "GEOObjects.h"

/**
 * \brief VtkConditionSource is a VTK source object for the visualization
 * of FEM conditions. As a vtkPolyDataAlgorithm it outputs polygonal data.
 */
class VtkConditionSource : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{

public:
	/// Create new objects with New() because of VTKs object reference counting.
	static VtkConditionSource* New();

	vtkTypeRevisionMacro(VtkConditionSource,vtkPolyDataAlgorithm);

	/// Sets the geo data as a vector
	void setData(const std::vector<GEOLIB::Point*>* points, const std::vector<GEOLIB::Polyline*>* lines, const std::vector<GEOLIB::Surface*>* surfaces,
		         std::vector<size_t> *pnt_idx, std::vector<size_t> *ply_idx, std::vector<size_t> *sfc_idx);

	/// Prints its data on a stream.
	void PrintSelf(ostream& os, vtkIndent indent);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	VtkConditionSource();
	~VtkConditionSource() {};

	/// Computes the polygonal data object.
	int RequestData(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

	int RequestInformation(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

	const std::vector<GEOLIB::Point*>* _points;
	const std::vector<GEOLIB::Polyline*>* _polylines;
	const std::vector<GEOLIB::Surface*>* _surfaces;
	std::vector<size_t>* _pnt_idx;
	std::vector<size_t>* _ply_idx;
	std::vector<size_t>* _sfc_idx;

private:

};

#endif // VTKCONDITIONSOURCE_H
