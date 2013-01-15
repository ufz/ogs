/**
 * \file   RandomEquationIdStorage.cpp
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Helper macros.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RandomEquationIdStorage.h"

namespace DiscreteLib
{

RandomEquationIdStorage::RandomEquationIdStorage(const std::vector<std::size_t> &sorted_list_pt_id) : _list_pt_id(sorted_list_pt_id)
{
    for (std::size_t i=0; i<sorted_list_pt_id.size(); i++)
        set(sorted_list_pt_id[i], 0);
}
RandomEquationIdStorage::RandomEquationIdStorage(std::size_t pt_id_start, std::size_t n)
{
    _list_pt_id.resize(n);
    for (std::size_t i=0; i<n; i++) {
        _list_pt_id[i] = pt_id_start + i;
        set(_list_pt_id[i], 0);
    }
}

bool RandomEquationIdStorage::hasKey(std::size_t pt_id) const
{
    return _list_pt_id.end() !=  std::find(_list_pt_id.begin(), _list_pt_id.end(), pt_id);
}

bool RandomEquationIdStorage::hasValue(std::size_t eqs_id) const
{
    return (_map_pt2eqs.right.find(eqs_id)!=_map_pt2eqs.right.end());
}


void RandomEquationIdStorage::key_range(std::size_t &i_start, std::size_t &i_end) const
{
    i_start = _list_pt_id[0];
    i_end = _list_pt_id.back()+1;
}

void RandomEquationIdStorage::activate(std::size_t pt_id, bool b)
{
    if (b) {
        _deactive.erase(pt_id);
    } else {
        _deactive.insert(pt_id);
        set(pt_id, -1);
    }
}

void RandomEquationIdStorage::set(std::size_t pt_id, long eqs_id)
{
    _map_pt2eqs.insert(value_type(pt_id, eqs_id));
}

std::size_t RandomEquationIdStorage::setAll(std::size_t dof_start, std::size_t delta_per_pt)
{
    std::size_t last = 0;
    if (_deactive.size()==0) {
        for (std::size_t i=0; i<_list_pt_id.size(); i++) {
            last = dof_start + i*delta_per_pt;
            set(_list_pt_id[i], last);
        }
    } else {
        std::size_t counter = 0;
        for (std::size_t i=0; i<_list_pt_id.size(); i++) {
            const std::size_t pt_id = _list_pt_id[i];
            if (_deactive.count(pt_id)>0) continue;
            last = dof_start + counter*delta_per_pt;
            set(pt_id, last);
            counter++;
        }
    }
    return last+1;
}

std::size_t RandomEquationIdStorage::address(std::size_t pt_id) const
{
    if (_map_pt2eqs.left.count(pt_id)==0) return -1;
    return _map_pt2eqs.left.find(pt_id)->second;
}

std::size_t RandomEquationIdStorage::key(std::size_t address_id) const
{
    if (_map_pt2eqs.right.count(address_id)==0) return -1;
    return _map_pt2eqs.right.find(address_id)->second;
}

} //end


