/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-26
 * \brief  Definition of the VtkMeshNodalCoordinatesTemplate class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "MeshLib/Node.h"

namespace MeshLib {

// Can't use vtkStandardNewMacro with a template.
template <class Scalar> VtkMeshNodalCoordinatesTemplate<Scalar> *
VtkMeshNodalCoordinatesTemplate<Scalar>::New()
{
    VTK_STANDARD_NEW_BODY(VtkMeshNodalCoordinatesTemplate<Scalar>)
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::PrintSelf(ostream &os, vtkIndent indent)
{
    this->VtkMeshNodalCoordinatesTemplate<Scalar>::Superclass::PrintSelf(
          os, indent);
    //os << indent << "XArray: " << this->XArray << std::endl;
    //os << indent << "YArray: " << this->YArray << std::endl;
    //os << indent << "ZArray: " << this->ZArray << std::endl;
    //os << indent << "TempDoubleArray: " << this->TempDoubleArray << std::endl;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::Initialize()
{
    this->_nodes = NULL;
    delete [] this->TempDoubleArray;
    this->TempDoubleArray = NULL;
    this->MaxId = -1;
    this->Size = 0;
    this->NumberOfComponents = 1;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::SetNodes(std::vector<MeshLib::Node*> const & nodes)
{
    Initialize();
    _nodes = &nodes;
    this->NumberOfComponents = 3;
    this->Size = this->NumberOfComponents * _nodes->size();
    this->MaxId = this->Size - 1;
    this->TempDoubleArray = new double [this->NumberOfComponents];
    this->Modified();
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::GetTuples(vtkIdList *ptIds, vtkAbstractArray *output)
{
    vtkDataArray *outArray = vtkDataArray::FastDownCast(output);
    if(!outArray)
    {
        vtkWarningMacro(<<"Input is not a vtkDataArray");
        return;
    }

    const vtkIdType numTuples = ptIds->GetNumberOfIds();

    outArray->SetNumberOfComponents(this->NumberOfComponents);
    outArray->SetNumberOfTuples(numTuples);

    const vtkIdType numPoints = ptIds->GetNumberOfIds();
    for(vtkIdType i = 0; i < numPoints; i++)
        outArray->SetTuple(i, this->GetTuple(ptIds->GetId(i)));
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::GetTuples(vtkIdType p1, vtkIdType p2, vtkAbstractArray *output)
{
    vtkDataArray *da = vtkDataArray::FastDownCast(output);
    if(!da)
    {
        vtkWarningMacro(<<"Input is not a vtkDataArray");
        return;
    }

    if(da->GetNumberOfComponents() != this->GetNumberOfComponents())
    {
        vtkErrorMacro(<<"Incorrect number of components in input array.");
        return;
    }

    for(vtkIdType daTubleId = 0; p1 <= p2; ++p1)
        da->SetTuple(daTubleId++, this->GetTuple(p1));
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::Squeeze()
{

}

template <class Scalar> vtkArrayIterator* VtkMeshNodalCoordinatesTemplate<Scalar>
::NewIterator()
{
    vtkErrorMacro(<<"Not implemented.");
    return NULL;
}

template <class Scalar> vtkIdType VtkMeshNodalCoordinatesTemplate<Scalar>
::LookupValue(vtkVariant value)
{
    bool valid = true;
    Scalar val = vtkVariantCast<Scalar>(value, &valid);
    if(valid)
        return this->Lookup(val, 0);
    return -1;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::LookupValue(vtkVariant value, vtkIdList *ids)
{
    bool valid = true;
    Scalar val = vtkVariantCast<Scalar>(value, &valid);
    ids->Reset();
    if(valid)
    {
        vtkIdType index = 0;
        while((index = this->Lookup(val, index)) >= 0)
            ids->InsertNextId(index++);
    }
}

template <class Scalar> vtkVariant VtkMeshNodalCoordinatesTemplate<Scalar>
::GetVariantValue(vtkIdType idx)
{
    return vtkVariant(this->GetValueReference(idx));
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::ClearLookup()
{
    // no fast lookup implemented
}

template <class Scalar> double* VtkMeshNodalCoordinatesTemplate<Scalar>
::GetTuple(vtkIdType i)
{
    this->GetTuple(i, this->TempDoubleArray);
    return this->TempDoubleArray;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::GetTuple(vtkIdType i, double *tuple)
{
    tuple[0] = (*(*this->_nodes)[i])[0];
    tuple[1] = (*(*this->_nodes)[i])[1];
    tuple[2] = (*(*this->_nodes)[i])[2];
}

template <class Scalar> vtkIdType VtkMeshNodalCoordinatesTemplate<Scalar>
::LookupTypedValue(Scalar value)
{
    return this->Lookup(value, 0);
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::LookupTypedValue(Scalar value, vtkIdList *ids)
{
    ids->Reset();
    vtkIdType index = 0;
    while((index = this->Lookup(value, index)) >= 0)
        ids->InsertNextId(index++);
}

template <class Scalar> Scalar& VtkMeshNodalCoordinatesTemplate<Scalar>
::GetValueReference(vtkIdType idx)
{
    const vtkIdType tuple = idx / this->NumberOfComponents;
    const vtkIdType comp = idx % this->NumberOfComponents;
    return (*(*this->_nodes)[tuple])[comp];
}

template <class Scalar> int VtkMeshNodalCoordinatesTemplate<Scalar>
::Allocate(vtkIdType, vtkIdType)
{
    vtkErrorMacro("Read only container.");
    return 0;
}

template <class Scalar> int VtkMeshNodalCoordinatesTemplate<Scalar>
::Resize(vtkIdType)
{
    vtkErrorMacro("Read only container.");
    return 0;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::SetNumberOfTuples(vtkIdType)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::SetTuple(vtkIdType, vtkIdType, vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::SetTuple(vtkIdType, const float *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::SetTuple(vtkIdType, const double *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertTuple(vtkIdType, vtkIdType, vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertTuple(vtkIdType, const float *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertTuple(vtkIdType, const double *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertTuples(vtkIdList *, vtkIdList *, vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertTuples(vtkIdType, vtkIdType, vtkIdType, vtkAbstractArray*)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> vtkIdType VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertNextTuple(vtkIdType, vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

template <class Scalar> vtkIdType VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertNextTuple(const float *)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

template <class Scalar> vtkIdType VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertNextTuple(const double *)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertVariantValue(vtkIdType /*idx*/, vtkVariant /*value*/)
{
    vtkErrorMacro("Read only container.");
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::DeepCopy(vtkAbstractArray *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::DeepCopy(vtkDataArray *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InterpolateTuple(vtkIdType, vtkIdList *, vtkAbstractArray *, double *)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InterpolateTuple(vtkIdType, vtkIdType, vtkAbstractArray*, vtkIdType,
                   vtkAbstractArray*, double)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::SetVariantValue(vtkIdType, vtkVariant)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::RemoveTuple(vtkIdType)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::RemoveFirstTuple()
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::RemoveLastTuple()
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::SetValue(vtkIdType, Scalar)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> vtkIdType VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertNextValue(Scalar)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertValue(vtkIdType, Scalar)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> VtkMeshNodalCoordinatesTemplate<Scalar>
::VtkMeshNodalCoordinatesTemplate()
    : _nodes(NULL), TempDoubleArray(NULL)
{

}

template <class Scalar> VtkMeshNodalCoordinatesTemplate<Scalar>
::~VtkMeshNodalCoordinatesTemplate()
{
    delete [] this->TempDoubleArray;
}

template <class Scalar> vtkIdType VtkMeshNodalCoordinatesTemplate<Scalar>
::Lookup(const Scalar &val, vtkIdType index)
{
    while(index <= this->MaxId)
    {
        if(this->GetValueReference(index++) == val)
            return index;
    }
    return -1;
}


#if !(VTK_MAJOR_VERSION < 7 || VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION < 1)
template <class Scalar> Scalar& VtkMeshNodalCoordinatesTemplate<Scalar>
::GetValueReference(vtkIdType idx) const
{
    const vtkIdType tuple = idx / this->NumberOfComponents;
    const vtkIdType comp = idx % this->NumberOfComponents;
    return (*(*this->_nodes)[tuple])[comp];
}

template <class Scalar> Scalar VtkMeshNodalCoordinatesTemplate<Scalar>
::GetValue(vtkIdType idx) const
{
    return this->GetValueReference(idx);
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::GetTypedTuple(vtkIdType tupleId, Scalar *tuple) const
{
    tuple[0] = (*(*this->_nodes)[tupleId])[0];
    tuple[1] = (*(*this->_nodes)[tupleId])[1];
    tuple[2] = (*(*this->_nodes)[tupleId])[2];
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::SetTypedTuple(vtkIdType, const Scalar*)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertTypedTuple(vtkIdType, const Scalar*)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> vtkIdType VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertNextTypedTuple(const Scalar *)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

#else
template <class Scalar> Scalar VtkMeshNodalCoordinatesTemplate<Scalar>
::GetValue(vtkIdType idx)
{
    return this->GetValueReference(idx);
}
template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::GetTupleValue(vtkIdType tupleId, Scalar *tuple)
{
    tuple[0] = (*(*this->_nodes)[tupleId])[0];
    tuple[1] = (*(*this->_nodes)[tupleId])[1];
    tuple[2] = (*(*this->_nodes)[tupleId])[2];
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::SetTupleValue(vtkIdType, const Scalar*)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> void VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertTupleValue(vtkIdType, const Scalar*)
{
    vtkErrorMacro("Read only container.");
    return;
}

template <class Scalar> vtkIdType VtkMeshNodalCoordinatesTemplate<Scalar>
::InsertNextTupleValue(const Scalar *)
{
    vtkErrorMacro("Read only container.");
    return -1;
}

#endif  // vtk version
} // end namespace
