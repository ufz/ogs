/**
 * \file
 * \author Lars Bilke
 * \date   2014-10-17
 * \brief  VtkMappedElementDataArrayTemplate is a adapter for cell data
 *         on elements of OGS meshes to VTK unstructured grids.
 *
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace InSituLib {

template <class Scalar>
class VtkMappedElementDataArrayTemplate:
	public vtkTypeTemplate<VtkMappedElementDataArrayTemplate<Scalar>,
	                       vtkMappedDataArray<Scalar> >
{
public:
	vtkMappedDataArrayNewInstanceMacro(VtkMappedElementDataArrayTemplate<Scalar>)
	static VtkMappedElementDataArrayTemplate *New();
	virtual void PrintSelf(std::ostream &os, vtkIndent indent);

	// Description:
	// Set the raw scalar arrays for the coordinate set. This class takes
	// ownership of the arrays and deletes them with delete[].
	void SetElements(std::vector<MeshLib::Element *> const * elements, vtkIdType numTuples);
	void SetElements(std::vector<MeshLib::Element *> const * elements, vtkIdType numTuples, bool save);

	// Reimplemented virtuals -- see superclasses for descriptions:
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

	// Description:
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
	vtkIdType InsertNextTuple(vtkIdType j, vtkAbstractArray *source);
	vtkIdType InsertNextTuple(const float *source);
	vtkIdType InsertNextTuple(const double *source);
	void DeepCopy(vtkAbstractArray *aa);
	void DeepCopy(vtkDataArray *da);
	void InterpolateTuple(vtkIdType i, vtkIdList *ptIndices,
	                      vtkAbstractArray* source, double* weights);
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
	VtkMappedElementDataArrayTemplate();
	~VtkMappedElementDataArrayTemplate();

	//Scalar *Array;
	std::vector<MeshLib::Element*> const * _elements;

private:
	VtkMappedElementDataArrayTemplate(
		const VtkMappedElementDataArrayTemplate &); // Not implemented.
	void operator=(
		const VtkMappedElementDataArrayTemplate &); // Not implemented.

	vtkIdType Lookup(const Scalar &val, vtkIdType startIndex);
	double TempDouble;
	// Description: If Save is true then this class won't delete that memory.
	// By default Save is false.
	bool Save;
};

} // end namespace InSituLib

#include "VtkMappedElementDataArrayTemplate-impl.h"

#endif // VTKMAPPEDELEMENTDATAARRAY_H_
