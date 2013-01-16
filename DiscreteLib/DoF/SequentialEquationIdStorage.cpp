/**
 * \file   SequentialEquationIdStorage.cpp
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

#include "SequentialEquationIdStorage.h"

#include <algorithm>
#include <cassert>

namespace DiscreteLib
{

bool SequentialEquationIdStorage::hasValue(std::size_t eqs_id) const
{
    bool in_range = (address(_pt_id_start) <= eqs_id && eqs_id <= address(_pt_id_start+_n-1));
    if (!in_range) return false;
    if ((eqs_id - _dof_begin)%_delta_per_pt!=0) return false;
    return true;
}

void SequentialEquationIdStorage::key_range(std::size_t &i_start, std::size_t &i_end) const
{
    i_start = _pt_id_start;
    i_end = i_start + _n;
}

void SequentialEquationIdStorage::activate(std::size_t pt_id, bool b)
{
    if (b) {
        _deactive.erase(pt_id);
    } else {
        _deactive.insert(pt_id);
    }
}

std::size_t SequentialEquationIdStorage::setAll(std::size_t address_start, std::size_t dn_pt)
{
    _dof_begin = address_start;
    _delta_per_pt = dn_pt;
    return _dof_begin + (_n-_deactive.size())*_delta_per_pt;
}

std::size_t SequentialEquationIdStorage::address(std::size_t pt_id) const
{
    assert(_pt_id_start<=pt_id && pt_id<_pt_id_start+_n);

    if (_deactive.count(pt_id)>0) return -1;

    if (_pt_id_start<=pt_id && pt_id<_pt_id_start+_n) {
        std::size_t loc = 0;
        if (_deactive.size()==0) {
            loc = _dof_begin + (pt_id-_pt_id_start)*_delta_per_pt;
        } else {
            std::size_t inactive_cnt = 0;
            for (std::size_t i=_pt_id_start; i<pt_id; i++) {
                if (_deactive.count(i)>0) inactive_cnt++;
            }
            loc = _dof_begin + (pt_id-_pt_id_start-inactive_cnt)*_delta_per_pt;
        }
        return loc;
    } else {
        return index_npos;
    }
}

} //end

