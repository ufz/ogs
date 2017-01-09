/**
 * \file
 * \author Lars Bilke
 * \date   2014-10-17
 * \brief  Implementation of the VtkMappedElementDataArrayTemplate class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <logog/include/logog.hpp>

#include <vtkIdList.h>
#include <vtkObjectFactory.h>
#include <vtkVariant.h>
#include <vtkVariantCast.h>

namespace MeshLib {

// Can't use vtkStandardNewMacro on a templated class.
template <class Scalar> VtkMappedPropertyVectorTemplate<Scalar> *
VtkMappedPropertyVectorTemplate<Scalar>::New()
{
    VTK_STANDARD_NEW_BODY(VtkMappedPropertyVectorTemplate<Scalar>)
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::PrintSelf(ostream &os, vtkIndent indent)
{
    this->VtkMappedPropertyVectorTemplate<Scalar>::Superclass::PrintSelf(
        os, indent);
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::SetPropertyVector(MeshLib::PropertyVector<Scalar> & propertyVector)
{
    this->Initialize();
    _propertyVector = &propertyVector;
    this->NumberOfComponents = _propertyVector->getNumberOfComponents();
    this->Size = this->NumberOfComponents *  _propertyVector->getNumberOfTuples();
    this->MaxId = this->Size - 1;
    this->Modified();
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::SetPropertyVector(MeshLib::PropertyVector<Scalar> & propertyVector, bool save)
{
    this->SetPropertyVector(propertyVector);
    this->Save = save;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::Initialize()
{
    this->MaxId = -1;
    this->Size = 0;
    this->NumberOfComponents = 1;
    // per default property vector deletion is managed elsewhere
    this->Save = true;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
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
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
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
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::Squeeze()
{
    // noop
}

//------------------------------------------------------------------------------
template <class Scalar> vtkArrayIterator*
VtkMappedPropertyVectorTemplate<Scalar>::NewIterator()
{
    vtkErrorMacro(<<"Not implemented.");
    return nullptr;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedPropertyVectorTemplate<Scalar>
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
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
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
template <class Scalar> vtkVariant VtkMappedPropertyVectorTemplate<Scalar>
::GetVariantValue(vtkIdType idx)
{
    return vtkVariant(this->GetValueReference(idx));
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::ClearLookup()
{
    // no-op, no fast lookup implemented.
}

//------------------------------------------------------------------------------
template <class Scalar> double* VtkMappedPropertyVectorTemplate<Scalar>
::GetTuple(vtkIdType i)
{
    this->TempDouble = (*this->_propertyVector)[i];
    return &this->TempDouble;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::GetTuple(vtkIdType i, double *tuple)
{
    *tuple = (*this->_propertyVector)[i];
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedPropertyVectorTemplate<Scalar>
::LookupTypedValue(Scalar value)
{
    return this->Lookup(value, 0);
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
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
template <class Scalar> Scalar& VtkMappedPropertyVectorTemplate<Scalar>
::GetValueReference(vtkIdType idx)
{
    // VTK has no concept of 'const', so we'll just cross our fingers
    // that no one writes to the returned reference.
    Scalar& value = const_cast<Scalar&>((*this->_propertyVector)[idx]);
    return value;
}

//------------------------------------------------------------------------------
template <class Scalar> int VtkMappedPropertyVectorTemplate<Scalar>
::Allocate(vtkIdType, vtkIdType)
{
    vtkErrorMacro("Read only container.");
    return 0;
}

//------------------------------------------------------------------------------
template <class Scalar> int VtkMappedPropertyVectorTemplate<Scalar>
::Resize(vtkIdType)
{
    vtkErrorMacro("Read only container.");
    return 0;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::SetNumberOfTuples(vtkIdType)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::SetTuple(vtkIdType, vtkIdType, vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::SetTuple(vtkIdType, const float *)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::SetTuple(vtkIdType, const double *)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::InsertTuple(vtkIdType, vtkIdType, vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::InsertTuple(vtkIdType, const float *)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::InsertTuple(vtkIdType, const double *)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::InsertTuples(vtkIdList *, vtkIdList *, vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
}

template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::InsertTuples(vtkIdType, vtkIdType, vtkIdType, vtkAbstractArray*)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedPropertyVectorTemplate<Scalar>
::InsertNextTuple(vtkIdType, vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedPropertyVectorTemplate<Scalar>
::InsertNextTuple(const float *)
{

    vtkErrorMacro("Read only container.");
    return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedPropertyVectorTemplate<Scalar>
::InsertNextTuple(const double *)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::DeepCopy(vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::DeepCopy(vtkDataArray *)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::InterpolateTuple(vtkIdType, vtkIdList *, vtkAbstractArray *, double *)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::InterpolateTuple(vtkIdType, vtkIdType, vtkAbstractArray*, vtkIdType,
                   vtkAbstractArray*, double)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::SetVariantValue(vtkIdType, vtkVariant)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::RemoveTuple(vtkIdType)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::RemoveFirstTuple()
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::RemoveLastTuple()
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::SetValue(vtkIdType, Scalar)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedPropertyVectorTemplate<Scalar>
::InsertNextValue(Scalar)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::InsertValue(vtkIdType, Scalar)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> void VtkMappedPropertyVectorTemplate<Scalar>
::InsertVariantValue(vtkIdType /*idx*/, vtkVariant /*value*/)
{
    vtkErrorMacro("Read only container.");
}

//------------------------------------------------------------------------------
template <class Scalar> VtkMappedPropertyVectorTemplate<Scalar>
::VtkMappedPropertyVectorTemplate()
  : Save(false)
{
}

//------------------------------------------------------------------------------
template <class Scalar> VtkMappedPropertyVectorTemplate<Scalar>
::~VtkMappedPropertyVectorTemplate()
{

}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType VtkMappedPropertyVectorTemplate<Scalar>
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

#if !(VTK_MAJOR_VERSION < 7 || VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION < 1)
//------------------------------------------------------------------------------
template <class Scalar>
Scalar VtkMappedPropertyVectorTemplate<Scalar>::GetValue(vtkIdType idx) const
{
    return (*this->_propertyVector)[idx];
}
//------------------------------------------------------------------------------
template <class Scalar>
void VtkMappedPropertyVectorTemplate<Scalar>::GetTypedTuple(vtkIdType tupleId,
                                                            Scalar* tuple) const
{
    *tuple = (*this->_propertyVector)[tupleId];
}
//------------------------------------------------------------------------------
template <class Scalar>
void VtkMappedPropertyVectorTemplate<Scalar>::SetTypedTuple(vtkIdType,
                                                            const Scalar*)
{
    vtkErrorMacro("Read only container.");
}
//------------------------------------------------------------------------------
template <class Scalar>
void VtkMappedPropertyVectorTemplate<Scalar>::InsertTypedTuple(vtkIdType,
                                                               const Scalar*)
{
    vtkErrorMacro("Read only container.");
}
//------------------------------------------------------------------------------
template <class Scalar>
vtkIdType VtkMappedPropertyVectorTemplate<Scalar>::InsertNextTypedTuple(
    const Scalar*)
{
    vtkErrorMacro("Read only container.");
    return -1;
}
#else
//------------------------------------------------------------------------------
template <class Scalar>
Scalar VtkMappedPropertyVectorTemplate<Scalar>::GetValue(vtkIdType idx)
{
    return (*this->_propertyVector)[idx];
}
//------------------------------------------------------------------------------
template <class Scalar>
void VtkMappedPropertyVectorTemplate<Scalar>::GetTupleValue(vtkIdType tupleId,
                                                            Scalar* tuple)
{
    *tuple = (*this->_propertyVector)[tupleId];
}
//------------------------------------------------------------------------------
template <class Scalar>
void VtkMappedPropertyVectorTemplate<Scalar>::SetTupleValue(vtkIdType,
                                                            const Scalar*)
{
    vtkErrorMacro("Read only container.");
}
//------------------------------------------------------------------------------
template <class Scalar>
void VtkMappedPropertyVectorTemplate<Scalar>::InsertTupleValue(vtkIdType,
                                                               const Scalar*)
{
    vtkErrorMacro("Read only container.");
}
//------------------------------------------------------------------------------
template <class Scalar>
vtkIdType VtkMappedPropertyVectorTemplate<Scalar>::InsertNextTupleValue(
    const Scalar*)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

#endif  // vtk version
} // end namespace MeshLib
