/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the ElementValueModification class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementValueModification.h"

#include <algorithm>

#include "logog/include/logog.hpp"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"


namespace MeshLib {

std::vector<unsigned> ElementValueModification::getMeshValues(const MeshLib::Mesh &mesh)
{
	const std::size_t nElements (mesh.getNElements());
	std::vector<unsigned> value_mapping;
	for (unsigned i=0; i<nElements; ++i)
	{
		bool exists(false);
		unsigned value (mesh.getElement(i)->getValue());
		const unsigned nValues (value_mapping.size());
		for (unsigned j=0; j<nValues; ++j)
		{
			if (value == value_mapping[j])
			{
				exists = true;
				break;
			}
		}
		if (!exists)
			value_mapping.push_back(value);
	}

	std::sort(value_mapping.begin(), value_mapping.end());
	return value_mapping;
}

bool ElementValueModification::replace(MeshLib::Mesh &mesh, unsigned old_value, unsigned new_value, bool replace_if_exists)
{
	std::vector<unsigned> value_mapping (ElementValueModification::getMeshValues(mesh));

	if (!replace_if_exists)
	{
		const unsigned nValues (value_mapping.size());
		for (unsigned j=0; j<nValues; ++j)
		{
			if (new_value == value_mapping[j])
			{
				WARN ("ElementValueModification::replaceElementValue() - Replacement value is already taken, no changes have been made.");
				return false;
			}
		}
	}
	const std::size_t nElements (mesh.getNElements());
	std::vector<MeshLib::Element*> &elements (const_cast<std::vector<MeshLib::Element*>&>(mesh.getElements()));
	for (unsigned i=0; i<nElements; ++i)
	{
		if (elements[i]->getValue() == old_value)
			elements[i]->setValue(new_value);
	}
	mesh.updateMaterialGroups();
	return true;
}

unsigned ElementValueModification::condense(MeshLib::Mesh &mesh)
{
	std::vector<unsigned> value_mapping (ElementValueModification::getMeshValues(mesh));
	std::vector<unsigned> reverse_mapping(value_mapping.back()+1, 0);
	const unsigned nValues (value_mapping.size());
	for (unsigned i=0; i<nValues; ++i)
		reverse_mapping[value_mapping[i]] = i;

	const std::size_t nElements (mesh.getNElements());
	std::vector<MeshLib::Element*> &elements (const_cast<std::vector<MeshLib::Element*>&>(mesh.getElements()));
	for (unsigned i=0; i<nElements; ++i)
		elements[i]->setValue(reverse_mapping[elements[i]->getValue()]);

	mesh.updateMaterialGroups();
	return nValues;
}

unsigned ElementValueModification::setByElementType(MeshLib::Mesh &mesh, MeshElemType ele_type, unsigned new_value)
{
	std::vector<MeshLib::Element*> &elements (const_cast<std::vector<MeshLib::Element*>&>(mesh.getElements()));
	unsigned nValues = 0;
	for (auto e : elements) {
		if (e->getGeomType()!=ele_type)
			continue;
		e->setValue(new_value);
		nValues++;
	}

	mesh.updateMaterialGroups();
	return nValues;
}

} // end namespace MeshLib
