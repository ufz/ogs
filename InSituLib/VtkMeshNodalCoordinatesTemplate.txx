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

#include "VtkMeshNodalCoordinatesTemplate.h"

#include <vtkIdList.h>
#include <vtkObjectFactory.h>
#include <vtkVariant.h>
#include <vtkVariantCast.h>

//------------------------------------------------------------------------------
// Can't use vtkStandardNewMacro with a template.
template <class Scalar> VtkMeshNodalCoordinatesTemplate<Scalar> *
VtkMeshNodalCoordinatesTemplate<Scalar>::New()
{
	VTK_STANDARD_NEW_BODY(VtkMeshNodalCoordinatesTemplate<Scalar>)
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::PrintSelf(ostream &os, vtkIndent indent)
{
	this->VtkMeshNodalCoordinatesTemplate<Scalar>::Superclass::PrintSelf(
	      os, indent);
	os << indent << "Blub" << std::endl;
	//os << indent << "XArray: " << this->XArray << std::endl;
	//os << indent << "YArray: " << this->YArray << std::endl;
	//os << indent << "ZArray: " << this->ZArray << std::endl;
	//os << indent << "TempDoubleArray: " << this->TempDoubleArray << std::endl;
}
