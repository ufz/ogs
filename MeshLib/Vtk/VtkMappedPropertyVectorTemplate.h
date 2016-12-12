/**
 * \file
 * \author Lars Bilke
 * \date   2014-10-17
 * \brief  VtkMappedElementDataArrayTemplate is a adapter for cell data
 *         on elements of OGS meshes to VTK unstructured grids.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKMAPPEDELEMENTDATAARRAY_H_
#define VTKMAPPEDELEMENTDATAARRAY_H_

#include <vtkMappedDataArray.h>
#include <vtkTypeTemplate.h>  // For templated vtkObject API
#include <vtkObjectFactory.h> // for vtkStandardNewMacro

#include "MeshLib/Elements/Element.h"

namespace MeshLib {

template <class Scalar>
class VtkMappedPropertyVectorTemplate :
    public vtkTypeTemplate<VtkMappedPropertyVectorTemplate<Scalar>,
                           vtkMappedDataArray<Scalar> >
{
public:
    vtkMappedDataArrayNewInstanceMacro(VtkMappedPropertyVectorTemplate<Scalar>);
    static VtkMappedPropertyVectorTemplate* New();
    virtual void PrintSelf(std::ostream &os, vtkIndent indent) override;

    // Description:
    // Set the raw scalar arrays for the coordinate set.
    void SetPropertyVector(MeshLib::PropertyVector<Scalar> & propertyVector);
    void SetPropertyVector(MeshLib::PropertyVector<Scalar> & propertyVector, bool save);

    // Reimplemented virtuals -- see superclasses for descriptions:
    void Initialize() override;
    void GetTuples(vtkIdList *ptIds, vtkAbstractArray *output) override;
    void GetTuples(vtkIdType p1, vtkIdType p2, vtkAbstractArray *output) override;
    void Squeeze() override;
    vtkArrayIterator *NewIterator() override;
    vtkIdType LookupValue(vtkVariant value) override;
    void LookupValue(vtkVariant value, vtkIdList *ids) override;
    vtkVariant GetVariantValue(vtkIdType idx) override;
    void ClearLookup() override;
    double* GetTuple(vtkIdType i) override;
    void GetTuple(vtkIdType i, double *tuple) override;
    vtkIdType LookupTypedValue(Scalar value) override;
    void LookupTypedValue(Scalar value, vtkIdList *ids) override;
    Scalar GetValue(vtkIdType idx);
    Scalar& GetValueReference(vtkIdType idx) override;
    void GetTupleValue(vtkIdType idx, Scalar *t);

    // Description:
    // This container is read only -- this method does nothing but print a
    // warning.
    int Allocate(vtkIdType sz, vtkIdType ext) override;
    int Resize(vtkIdType numTuples) override;
    void SetNumberOfTuples(vtkIdType number) override;
    void SetTuple(vtkIdType i, vtkIdType j, vtkAbstractArray *source) override;
    void SetTuple(vtkIdType i, const float *source) override;
    void SetTuple(vtkIdType i, const double *source) override;
    void InsertTuple(vtkIdType i, vtkIdType j, vtkAbstractArray *source) override;
    void InsertTuple(vtkIdType i, const float *source) override;
    void InsertTuple(vtkIdType i, const double *source) override;
    void InsertTuples(vtkIdList *dstIds, vtkIdList *srcIds,
                      vtkAbstractArray *source) override;
    void InsertTuples(vtkIdType, vtkIdType, vtkIdType, vtkAbstractArray*) override;
    vtkIdType InsertNextTuple(vtkIdType j, vtkAbstractArray *source) override;
    vtkIdType InsertNextTuple(const float *source) override;
    vtkIdType InsertNextTuple(const double *source) override;
    void InsertVariantValue(vtkIdType idx, vtkVariant value) override;
    void DeepCopy(vtkAbstractArray *aa) override;
    void DeepCopy(vtkDataArray *da) override;
    void InterpolateTuple(vtkIdType i, vtkIdList *ptIndices,
                          vtkAbstractArray* source, double* weights) override;
    void InterpolateTuple(vtkIdType i, vtkIdType id1, vtkAbstractArray *source1,
                          vtkIdType id2, vtkAbstractArray *source2, double t) override;
    void SetVariantValue(vtkIdType idx, vtkVariant value) override;
    void RemoveTuple(vtkIdType id) override;
    void RemoveFirstTuple() override;
    void RemoveLastTuple() override;
    void SetTupleValue(vtkIdType i, const Scalar *t);
    void InsertTupleValue(vtkIdType i, const Scalar *t);
    vtkIdType InsertNextTupleValue(const Scalar *t);
    void SetValue(vtkIdType idx, Scalar value) override;
    vtkIdType InsertNextValue(Scalar v) override;
    void InsertValue(vtkIdType idx, Scalar v) override;

protected:
    VtkMappedPropertyVectorTemplate();
    ~VtkMappedPropertyVectorTemplate();

    MeshLib::PropertyVector<Scalar> * _propertyVector;

private:
    VtkMappedPropertyVectorTemplate(
        const VtkMappedPropertyVectorTemplate &); // Not implemented.
    void operator=(
        const VtkMappedPropertyVectorTemplate &); // Not implemented.

    vtkIdType Lookup(const Scalar &val, vtkIdType startIndex);
    double TempDouble;
    // Description: If Save is true then this class won't delete that memory.
    // By default Save is false.
    bool Save;
};

} // end namespace MeshLib

#include "VtkMappedPropertyVectorTemplate-impl.h"

#endif // VTKMAPPEDELEMENTDATAARRAY_H_
