/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-26
 * \brief  Definition of the ElementStatus class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKMESHNODALCOORDINATES_H_
#define VTKMESHNODALCOORDINATES_H_

#include <vtkMappedDataArray.h>
#include <vtkTypeTemplate.h>    // For templated vtkObject API
#include <vtkObjectFactory.h>   // for vtkStandardNewMacro

template <class Scalar>
class VtkMeshNodalCoordinatesTemplate:
	public vtkTypeTemplate<VtkMeshNodalCoordinatesTemplate<Scalar>,
	                       vtkMappedDataArray<Scalar> >
{
public:
	vtkMappedDataArrayNewInstanceMacro(
		VtkMeshNodalCoordinatesTemplate<Scalar>)
	static VtkMeshNodalCoordinatesTemplate *New();
	virtual void PrintSelf(ostream &os, vtkIndent indent);

	void SetNodes(std::vector<MeshLib::Node*> const & nodes);
};

#include "VtkMeshNodalCoordinatesTemplate.txx"

#endif // VTKMESHNODALCOORDINATES_H_
