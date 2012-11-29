/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkPointsSource.h
 *
 * Created on 2010-02-03 by Lars Bilke
 */

#ifndef VTKPOINTSSOURCE_H
#define VTKPOINTSSOURCE_H

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include <vtkPolyDataAlgorithm.h>

#include "GEOObjects.h"

/**
 * \brief VtkPointsSource is a VTK source object for the visualization
 * of point data. As a vtkPolyDataAlgorithm it outputs polygonal data.
 */
class VtkPointsSource : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
	/// Create new objects with New() because of VTKs object reference counting.
	static VtkPointsSource* New();

	vtkTypeRevisionMacro(VtkPointsSource,vtkPolyDataAlgorithm);

	/// Sets the points as a vector
	void setPoints(const std::vector<GeoLib::Point*>* points) { _points = points; }

	/// Prints its data on a stream.
	void PrintSelf(ostream& os, vtkIndent indent);

	virtual void SetUserProperty(QString name, QVariant value);

protected:
	VtkPointsSource();
	~VtkPointsSource() {}

	/// Computes the polygonal data object.
	int RequestData(vtkInformation* request,
	                vtkInformationVector** inputVector,
	                vtkInformationVector* outputVector);

	int RequestInformation(vtkInformation* request,
	                       vtkInformationVector** inputVector,
	                       vtkInformationVector* outputVector);

	/// The points to visualize
	const std::vector<GeoLib::Point*>* _points;

private:
};

#endif // VTKPOINTSSOURCE_H
