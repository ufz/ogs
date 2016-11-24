/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-26
 * \brief  VtkMeshNodalCoordinatesTemplate is a adapter for node coordinates of
 *         OGS meshes to VTK unstructured grids.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace MeshLib
{
    class Node;
}

namespace MeshLib
{

template <class Scalar>
class VtkMeshNodalCoordinatesTemplate:
    public vtkTypeTemplate<VtkMeshNodalCoordinatesTemplate<Scalar>,
                           vtkMappedDataArray<Scalar> >
{
public:
    vtkMappedDataArrayNewInstanceMacro(VtkMeshNodalCoordinatesTemplate<Scalar>);
    static VtkMeshNodalCoordinatesTemplate *New();
    virtual void PrintSelf(std::ostream &os, vtkIndent indent);

    /// Pass the nodes from OGS mesh
    void SetNodes(std::vector<MeshLib::Node*> const & nodes);

    // Reimplemented virtuals -- see superclasses for descriptions
    void Initialize();
    void GetTuples(vtkIdList *ptIds, vtkAbstractArray *output);
    void GetTuples(vtkIdType p1, vtkIdType p2, vtkAbstractArray *output);
    void Squeeze();
    vtkArrayIterator *NewIterator();
    vtkIdType LookupValue(vtkVariant value);
    void LookupValue(vtkVariant value, vtkIdList *ids);
    vtkVariant GetVariantValue(vtkIdType idx);
    void ClearLookup();
    double* GetTuple(vtkIdType i);
    void GetTuple(vtkIdType i, double *tuple);
    vtkIdType LookupTypedValue(Scalar value);
    void LookupTypedValue(Scalar value, vtkIdList *ids);
    Scalar GetValue(vtkIdType idx);
    Scalar& GetValueReference(vtkIdType idx);
    void GetTupleValue(vtkIdType idx, Scalar *t);

    // This container is read only -- this method does nothing but print a
    // warning.
    int Allocate(vtkIdType sz, vtkIdType ext);
    int Resize(vtkIdType numTuples);
    void SetNumberOfTuples(vtkIdType number);
    void SetTuple(vtkIdType i, vtkIdType j, vtkAbstractArray *source);
    void SetTuple(vtkIdType i, const float *source);
    void SetTuple(vtkIdType i, const double *source);
    void InsertTuple(vtkIdType i, vtkIdType j, vtkAbstractArray *source);
    void InsertTuple(vtkIdType i, const float *source);
    void InsertTuple(vtkIdType i, const double *source);
    void InsertTuples(vtkIdList *dstIds, vtkIdList *srcIds,
                      vtkAbstractArray *source);
    void InsertTuples(vtkIdType, vtkIdType, vtkIdType, vtkAbstractArray*);
    vtkIdType InsertNextTuple(vtkIdType j, vtkAbstractArray *source);
    vtkIdType InsertNextTuple(const float *source);
    vtkIdType InsertNextTuple(const double *source);
    void InsertVariantValue(vtkIdType idx, vtkVariant value);
    void DeepCopy(vtkAbstractArray *aa);
    void DeepCopy(vtkDataArray *da);
    void InterpolateTuple(vtkIdType i, vtkIdList *ptIndices,
                          vtkAbstractArray* source,  double* weights);
    void InterpolateTuple(vtkIdType i, vtkIdType id1, vtkAbstractArray *source1,
                          vtkIdType id2, vtkAbstractArray *source2, double t);
    void SetVariantValue(vtkIdType idx, vtkVariant value);
    void RemoveTuple(vtkIdType id);
    void RemoveFirstTuple();
    void RemoveLastTuple();
    void SetTupleValue(vtkIdType i, const Scalar *t);
    void InsertTupleValue(vtkIdType i, const Scalar *t);
    vtkIdType InsertNextTupleValue(const Scalar *t);
    void SetValue(vtkIdType idx, Scalar value);
    vtkIdType InsertNextValue(Scalar v);
    void InsertValue(vtkIdType idx, Scalar v);


protected:
    VtkMeshNodalCoordinatesTemplate();
    ~VtkMeshNodalCoordinatesTemplate();

    const std::vector<MeshLib::Node*>* _nodes;

private:
    // Not implemented
    VtkMeshNodalCoordinatesTemplate(const VtkMeshNodalCoordinatesTemplate &);
    void operator=(const VtkMeshNodalCoordinatesTemplate &);

    vtkIdType Lookup(const Scalar &val, vtkIdType startIndex);
    double *TempDoubleArray;
};

} // end namespace

#include "VtkMeshNodalCoordinatesTemplate-impl.h"

#endif // VTKMESHNODALCOORDINATES_H_
