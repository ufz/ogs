/**
 * \file VtkConditionSource.h
 * 2011/03/02 KR Initial implementation
 *
 */

#ifndef VTKCONDITIONSOURCE_H
#define VTKCONDITIONSOURCE_H

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include <vtkPolyDataAlgorithm.h>

#include "GEOObjects.h"
//#include <vtkSmartPointer.h>

class FEMCondition;
//class vtkPoints;
//class vtkDoubleArray;

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

	/// Sets the FEMCondition that need to be visualised. The geometry points array is needed because polylines and surfaces are linked to this data.
	void setData(const std::vector<GeoLib::Point*>* points,
	             const std::vector<FEMCondition*>* conds);

	/// Prints its data on a stream.
	void PrintSelf(ostream& os, vtkIndent indent);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	VtkConditionSource();
	~VtkConditionSource();

	/// Computes the polygonal data object.
	int RequestData(vtkInformation* request,
	                vtkInformationVector** inputVector,
	                vtkInformationVector* outputVector);

	int RequestInformation(vtkInformation* request,
	                       vtkInformationVector** inputVector,
	                       vtkInformationVector* outputVector);

private:
	//size_t getIndex(size_t idx, vtkSmartPointer<vtkPoints> newPoints, vtkSmartPointer<vtkDoubleArray> scalars, std::map<size_t, size_t> &idx_map);

	const std::vector<GeoLib::Point*>* _points;
	const std::vector<FEMCondition*>* _cond_vec;
	std::map<FiniteElement::DistributionType, vtkIdType> _dis_type_map;
};

#endif // VTKCONDITIONSOURCE_H
