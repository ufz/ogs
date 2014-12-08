/**
 * \file
 * \brief  Definition of the class Properties that implements a container of
 *         properties.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROPERTYVECTOR_H_
#define PROPERTYVECTOR_H_

#include <algorithm>
#include <vector>
#include <type_traits>
#include <memory>

namespace MeshLib
{

template <typename PROP_VAL_TYPE>
class PropertyVector : public std::vector<PROP_VAL_TYPE>
{
public:
	PropertyVector()
		: std::vector<PROP_VAL_TYPE>()
	{}

	explicit PropertyVector(std::size_t size)
		: std::vector<PROP_VAL_TYPE>(size)
	{}
};

template <typename T>
class PropertyVector<T*> : public std::vector<T*>
{
public:
	PropertyVector(std::size_t n_mat_groups,
		std::vector<std::size_t> const& mat_group_idx_map)
		: std::vector<T*>(n_mat_groups), _mat_group_idx_map(mat_group_idx_map)
	{}

	~PropertyVector()
	{
		std::for_each(
			this->begin(), this->end(), std::default_delete<T>()
		);
	}

	T* const& operator[](std::size_t id) const
	{
		return (*static_cast<std::vector<T*> const *>(this))[_mat_group_idx_map[id]];
	}

private:
	std::vector<std::size_t> _mat_group_idx_map;
};

} // end namespace MeshLib

#endif
