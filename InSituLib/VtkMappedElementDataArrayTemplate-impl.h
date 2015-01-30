/**
 * \file
 * \author Lars Bilke
 * \date   2014-10-17
 * \brief  Implementation of the VtkMappedElementDataArrayTemplate class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "logog/include/logog.hpp"

#include <vtkIdList.h>
#include <vtkObjectFactory.h>
#include <vtkVariant.h>
#include <vtkVariantCast.h>

namespace InSituLib {

// Can't use vtkStandardNewMacro on a templated class.
template <class Scalar> VtkMappedElementDataArrayTemplate<Scalar> *
VtkMappedElementDataArrayTemplate<Scalar>::New()
{
	VTK_STANDARD_NEW_BODY(VtkMappedElementDataArrayTemplate<Scalar>)
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::PrintSelf(ostream &os, vtkIndent indent)
{
	this->VtkMappedElementDataArrayTemplate<Scalar>::Superclass::PrintSelf(
		os, indent);
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::SetPropertyVector(MeshLib::PropertyVector<Scalar> & propertyVector)
{
	this->Initialize();
	_propertyVector = &propertyVector;
	this->NumberOfComponents = _propertyVector->getTupleSize();
	this->Size = this->NumberOfComponents *  _propertyVector->size();
	this->MaxId = this->Size - 1;
	this->Modified();
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::SetPropertyVector(MeshLib::PropertyVector<Scalar> & propertyVector, bool save)
{
	this->SetPropertyVector(propertyVector);
	this->Save = save;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::Initialize()
{
	this->MaxId = -1;
	this->Size = 0;
	this->NumberOfComponents = 1;
	// per default property vector deletion is managed elsewhere
	this->Save = true;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::GetTuples(vtkIdList *ptIds, vtkAbstractArray *output)
{
	vtkDataArray *da = vtkDataArray::FastDownCast(output);
	if (!da)
	{
		WARN("VtkMappedElementDataArrayTemplate::GetTuples(): Input is not a vtkDataArray");
		return;
	}

	if (da->GetNumberOfComponents() != this->GetNumberOfComponents())
	{
		WARN("VtkMappedElementDataArrayTemplate::GetTuples(): Incorrect number of components in input array.");
		return;
	}

	const vtkIdType numPoints = ptIds->GetNumberOfIds();
	for (vtkIdType i = 0; i < numPoints; ++i)
	{
		da->SetTuple(i, this->GetTuple(ptIds->GetId(i)));
	}
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::GetTuples(vtkIdType p1, vtkIdType p2, vtkAbstractArray *output)
{
	vtkDataArray *da = vtkDataArray::FastDownCast(output);
	if (!da)
	{
		vtkErrorMacro(<<"Input is not a vtkDataArray");
		return;
	}

	if (da->GetNumberOfComponents() != this->GetNumberOfComponents())
	{
		vtkErrorMacro(<<"Incorrect number of components in input array.");
		return;
	}

	for (vtkIdType daTupleId = 0; p1 <= p2; ++p1)
	{
		da->SetTuple(daTupleId++, this->GetTuple(p1));
	}
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::Squeeze()
{
	// noop
}

//------------------------------------------------------------------------------
template <class Scalar> vtkArrayIterator*
VtkMappedElementDataArrayTemplate<Scalar>::NewIterator()
{
	vtkErrorMacro(<<"Not implemented.");
	return nullptr;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedElementDataArrayTemplate<Scalar>
::LookupValue(vtkVariant value)
{
	bool valid = true;
	Scalar val = vtkVariantCast<Scalar>(value, &valid);
	if (valid)
	{
		return this->Lookup(val, 0);
	}
	return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::LookupValue(vtkVariant value, vtkIdList *ids)
{
	bool valid = true;
	Scalar val = vtkVariantCast<Scalar>(value, &valid);
	ids->Reset();
	if (valid)
	{
		vtkIdType index = 0;
		while ((index = this->Lookup(val, index)) >= 0)
		{
			ids->InsertNextId(index);
			++index;
		}
	}
}

//------------------------------------------------------------------------------
template <class Scalar> vtkVariant VtkMappedElementDataArrayTemplate<Scalar>
::GetVariantValue(vtkIdType idx)
{
	return vtkVariant(this->GetValueReference(idx));
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::ClearLookup()
{
	// no-op, no fast lookup implemented.
}

//------------------------------------------------------------------------------
template <class Scalar> double* VtkMappedElementDataArrayTemplate<Scalar>
::GetTuple(vtkIdType i)
{
	this->TempDouble = (*this->_propertyVector)[i];
	return &this->TempDouble;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::GetTuple(vtkIdType i, double *tuple)
{
	*tuple = (*this->_propertyVector)[i];
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedElementDataArrayTemplate<Scalar>
::LookupTypedValue(Scalar value)
{
	return this->Lookup(value, 0);
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::LookupTypedValue(Scalar value, vtkIdList *ids)
{
	ids->Reset();
	vtkIdType index = 0;
	while ((index = this->Lookup(value, index)) >= 0)
	{
		ids->InsertNextId(index);
		++index;
	}
}

//------------------------------------------------------------------------------
template <class Scalar> Scalar VtkMappedElementDataArrayTemplate<Scalar>
::GetValue(vtkIdType idx)
{
	return (*this->_propertyVector)[idx];
}

//------------------------------------------------------------------------------
template <class Scalar> Scalar& VtkMappedElementDataArrayTemplate<Scalar>
::GetValueReference(vtkIdType idx)
{
	// VTK has no concept of 'const', so we'll just cross our fingers
	// that no one writes to the returned reference.
	Scalar& value = const_cast<Scalar&>((*this->_propertyVector)[idx]);
	return value;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::GetTupleValue(vtkIdType tupleId, Scalar *tuple)
{
	*tuple = (*this->_propertyVector)[tupleId];
}

//------------------------------------------------------------------------------
template <class Scalar> int VtkMappedElementDataArrayTemplate<Scalar>
::Allocate(vtkIdType, vtkIdType)
{
	vtkErrorMacro("Read only container.")
	return 0;
}

//------------------------------------------------------------------------------
template <class Scalar> int VtkMappedElementDataArrayTemplate<Scalar>
::Resize(vtkIdType)
{
	vtkErrorMacro("Read only container.")
	return 0;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::SetNumberOfTuples(vtkIdType)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::SetTuple(vtkIdType, vtkIdType, vtkAbstractArray *)
{
	vtkErrorMacro("Read only container.")}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::SetTuple(vtkIdType, const float *)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::SetTuple(vtkIdType, const double *)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::InsertTuple(vtkIdType, vtkIdType, vtkAbstractArray *)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::InsertTuple(vtkIdType, const float *)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::InsertTuple(vtkIdType, const double *)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::InsertTuples(vtkIdList *, vtkIdList *, vtkAbstractArray *)
{
	vtkErrorMacro("Read only container.")
}

template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::InsertTuples(vtkIdType, vtkIdType, vtkIdType, vtkAbstractArray*)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedElementDataArrayTemplate<Scalar>
::InsertNextTuple(vtkIdType, vtkAbstractArray *)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedElementDataArrayTemplate<Scalar>
::InsertNextTuple(const float *)
{

	vtkErrorMacro("Read only container.")
	return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedElementDataArrayTemplate<Scalar>
::InsertNextTuple(const double *)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::DeepCopy(vtkAbstractArray *)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::DeepCopy(vtkDataArray *)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::InterpolateTuple(vtkIdType, vtkIdList *, vtkAbstractArray *, double *)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::InterpolateTuple(vtkIdType, vtkIdType, vtkAbstractArray*, vtkIdType,
				   vtkAbstractArray*, double)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::SetVariantValue(vtkIdType, vtkVariant)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::RemoveTuple(vtkIdType)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::RemoveFirstTuple()
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::RemoveLastTuple()
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::SetTupleValue(vtkIdType, const Scalar*)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::InsertTupleValue(vtkIdType, const Scalar*)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedElementDataArrayTemplate<Scalar>
::InsertNextTupleValue(const Scalar *)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::SetValue(vtkIdType, Scalar)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedElementDataArrayTemplate<Scalar>
::InsertNextValue(Scalar)
{
	vtkErrorMacro("Read only container.")
	return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedElementDataArrayTemplate<Scalar>
::InsertValue(vtkIdType, Scalar)
{
	vtkErrorMacro("Read only container.")
}

//------------------------------------------------------------------------------
template <class Scalar> VtkMappedElementDataArrayTemplate<Scalar>
::VtkMappedElementDataArrayTemplate()
  : Save(false)
{
}

//------------------------------------------------------------------------------
template <class Scalar> VtkMappedElementDataArrayTemplate<Scalar>
::~VtkMappedElementDataArrayTemplate()
{

}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedElementDataArrayTemplate<Scalar>
::Lookup(const Scalar &val, vtkIdType index)
{
	while (index <= this->MaxId)
	{
		if (this->GetValueReference(index++) == val)
		{
			return index;
		}
	}
	return -1;
}

} // end namespace InSituLib
