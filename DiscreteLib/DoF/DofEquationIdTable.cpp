/**
 * \file   DofEquationIdTable.cpp
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
#include "DofEquationIdTable.h"
#include <limits>

namespace DiscreteLib
{

///
DofEquationIdTable::~DofEquationIdTable()
{
    for (std::size_t i=0; i<_map_var2dof.size(); i++) {
        BaseLib::releaseObjectsInStdMap(_map_var2dof[i]);
    }
    _map_var2dof.clear();
}

std::size_t DofEquationIdTable::addVariableDoFs(std::size_t mesh_id, const std::vector<std::size_t> &list_dof_pt_id)
{
    RandomEquationIdStorage* address = new RandomEquationIdStorage(list_dof_pt_id);
    std::size_t var_id = _map_var2dof.size();
    _map_var2dof.resize(_map_var2dof.size()+1);
    _map_var2dof[var_id][mesh_id] = address;
    return var_id;
}

std::size_t DofEquationIdTable::addVariableDoFs(std::size_t mesh_id, std::size_t dof_pt_id_begin, std::size_t dof_pt_count)
{
    RandomEquationIdStorage* address = new RandomEquationIdStorage(dof_pt_id_begin, dof_pt_count);
    std::size_t var_id = _map_var2dof.size();
    _map_var2dof.resize(_map_var2dof.size()+1);
    _map_var2dof[var_id][mesh_id] = address;
    _map_msh2var[mesh_id].push_back(var_id);
    return var_id;
}

void DofEquationIdTable::addVariableDoFs(std::size_t n_var, std::size_t mesh_id, std::size_t dof_pt_id_begin, std::size_t dof_pt_count)
{
    _is_seq_address = true;
    for (std::size_t i=0; i<n_var; i++) {
        SequentialEquationIdStorage* address = new SequentialEquationIdStorage(dof_pt_id_begin, dof_pt_count);
        std::size_t var_id = _map_var2dof.size();
        _map_var2dof.resize(_map_var2dof.size()+1);
        _map_var2dof[var_id][mesh_id] = address;
        _map_msh2var[mesh_id].push_back(var_id);
    }
}

void DofEquationIdTable::deactivateDoFs(std::size_t var_id, std::size_t mesh_id, std::size_t pt_id)
{
    getPointEquationIdTable(var_id, mesh_id)->activate(pt_id, false);
}

void DofEquationIdTable::setGhostPoints(std::size_t mesh_id, std::vector<std::size_t> &list_pt_id)
{
    _ghost_pt[mesh_id].assign(list_pt_id.begin(), list_pt_id.end());
}


bool DofEquationIdTable::isGhostPoint(std::size_t mesh_id, std::size_t pt_id) const
{
    if (_ghost_pt.size()==0) return false;
    assert(_ghost_pt.count(mesh_id)>0);
    const std::vector<std::size_t>& ghosts = _ghost_pt.find(mesh_id)->second;
    return std::find(ghosts.begin(), ghosts.end(), pt_id)!=ghosts.end();
}

std::size_t DofEquationIdTable::getTotalNumberOfGhostPoints() const
{
    std::size_t count = 0;
    for (std::map<std::size_t, std::vector<std::size_t> >::const_iterator itr=_ghost_pt.begin(); itr!=_ghost_pt.end(); ++itr) {
        count += itr->second.size();
    }
    return count;
}

bool DofEquationIdTable::isActiveDoF(std::size_t var_id, std::size_t mesh_id, std::size_t pt_id) const
{
    return getPointEquationIdTable(var_id, mesh_id)->isActive(pt_id);
}

std::size_t DofEquationIdTable::mapEqsID(std::size_t var_id, std::size_t mesh_id, std::size_t pt_id) const
{
    const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
    return add->address(pt_id);
}

void DofEquationIdTable::mapEqsID(std::size_t var_id, std::size_t mesh_id, const std::vector<std::size_t> &pt_id, std::vector<std::size_t> &eqs_id) const
{
    eqs_id.resize(pt_id.size());
    const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
    for (std::size_t i=0; i<pt_id.size(); i++)
        eqs_id[i] = add->address(pt_id[i]);
}

void DofEquationIdTable::mapEqsID(std::size_t mesh_id, const std::vector<std::size_t> &pt_id, std::vector<std::size_t> &eqs_id) const
{
    if (_map_msh2var.count(mesh_id)==0) return;

    const std::vector<std::size_t> &list_var = _map_msh2var.find(mesh_id)->second;
    const std::size_t n_pt = pt_id.size();
    const std::size_t n_dof_per_pt = list_var.size();
    const std::size_t n_total_dof = n_pt * n_dof_per_pt;
    eqs_id.resize(n_total_dof);

    if (_local_numbering_type==DofNumberingType::BY_VARIABLE) {
        for (std::size_t i=0; i<list_var.size(); i++) {
            std::size_t var_id = list_var[i];
            const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
            for (std::size_t j=0; j<pt_id.size(); j++) {
                eqs_id[j+i*n_pt] = add->address(pt_id[j]);
            }
        }
    } else {
        for (std::size_t i=0; i<list_var.size(); i++) {
            std::size_t var_id = list_var[i];
            const IEquationIdStorage* add = getPointEquationIdTable(var_id, mesh_id);
            for (std::size_t j=0; j<pt_id.size(); j++) {
                eqs_id[i + j*n_dof_per_pt] = add->address(pt_id[j]);
            }
        }
    }
}

void DofEquationIdTable::mapEqsID(
        std::size_t mesh_id, const std::vector<std::size_t> &list_pt_id,
        std::vector<std::size_t> &list_eqs_id,
        std::vector<std::size_t> &list_eqs_id_without_ghost
        ) const
{
    if (_map_msh2var.count(mesh_id)==0) return;

    const std::vector<std::size_t> &list_var = _map_msh2var.find(mesh_id)->second;
    const std::size_t n_pt = list_pt_id.size();
    const std::size_t n_dof_per_pt = list_var.size();
    const std::size_t n_total_dof = n_pt * n_dof_per_pt;
    list_eqs_id.resize(n_total_dof);
    list_eqs_id_without_ghost.resize(n_total_dof);

    for (std::size_t i_var=0; i_var<list_var.size(); i_var++) {
        const std::size_t var_id = list_var[i_var];
        const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
        for (std::size_t i_pt=0; i_pt<n_pt; i_pt++) {
            const std::size_t pt_id = list_pt_id[i_pt];
            std::size_t pos;
            if (_local_numbering_type==DofNumberingType::BY_VARIABLE) {
                pos = i_pt + i_var * n_pt;
            } else {
                pos = i_var + i_pt * n_dof_per_pt;
            }
            std::size_t this_eqs_id = table->address(pt_id);
            list_eqs_id[pos] = isGhostPoint(mesh_id, pt_id) ? IEquationIdStorage::index_npos : this_eqs_id;
            list_eqs_id_without_ghost[pos] = this_eqs_id;
        }
    }
}

void DofEquationIdTable::mapEqsIDreduced(
        std::size_t mesh_id, const std::vector<std::size_t> &list_pt_id,
        std::vector<std::size_t> &list_eqs_id,
        std::vector<std::size_t> &list_eqs_id_without_ghost
        ) const
{
    if (_map_msh2var.count(mesh_id)==0) return;

    const std::vector<std::size_t> &list_var = _map_msh2var.find(mesh_id)->second;
    const std::size_t n_pt = list_pt_id.size();
    const std::size_t n_dof_per_pt = list_var.size();
    const std::size_t n_total_dof = n_pt * n_dof_per_pt;
    list_eqs_id.reserve(n_total_dof);
    list_eqs_id_without_ghost.reserve(n_total_dof);

    if (_local_numbering_type==DofNumberingType::BY_VARIABLE) {
        for (std::size_t i_var=0; i_var<list_var.size(); i_var++) {
            const std::size_t var_id = list_var[i_var];
            const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
            std::vector<std::size_t> active_pt_list;
            std::vector<std::size_t> active_eqs_id_list;
            for (std::size_t i_pt=0; i_pt<n_pt; i_pt++) {
                const std::size_t pt_id = list_pt_id[i_pt];
                if (!table->hasKey(pt_id)) continue;
                std::size_t this_eqs_id = table->address(pt_id);
                list_eqs_id.push_back(isGhostPoint(mesh_id, pt_id) ? IEquationIdStorage::index_npos : this_eqs_id);
                list_eqs_id_without_ghost.push_back(this_eqs_id);
            }
        }
    } else { // by point
        std::vector<std::vector<std::size_t> > var_active_pt_list(list_var.size());
        std::vector<std::vector<std::size_t> > var_active_eqs_id_list(list_var.size());
        for (std::size_t i_pt=0; i_pt<n_pt; i_pt++) {
            const std::size_t pt_id = list_pt_id[i_pt];
            for (std::size_t i_var=0; i_var<list_var.size(); i_var++) {
                const std::size_t var_id = list_var[i_var];
                const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
                if (!table->hasKey(pt_id)) continue;
                std::size_t this_eqs_id = table->address(pt_id);
                list_eqs_id.push_back(isGhostPoint(mesh_id, pt_id) ? IEquationIdStorage::index_npos : this_eqs_id);
                list_eqs_id_without_ghost.push_back(this_eqs_id);
            }
        }
    }
}

void DofEquationIdTable::createLocalMappingTable(
    std::size_t mesh_id, const std::vector<std::size_t> &list_pt_id,
    DofEquationIdTable &local_table
    ) const
{
    const std::vector<std::size_t> &list_var = _map_msh2var.find(mesh_id)->second;
    const std::size_t n_pt = list_pt_id.size();

    if (_local_numbering_type==DofNumberingType::BY_VARIABLE) {
        std::size_t active_entry_counter = 0;
        for (std::size_t i_var=0; i_var<list_var.size(); i_var++) {
            const std::size_t var_id = list_var[i_var];
            const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
            std::vector<std::size_t> active_pt_list;
            std::vector<std::size_t> active_eqs_id_list;
            for (std::size_t i_pt=0; i_pt<n_pt; i_pt++) {
                const std::size_t pt_id = list_pt_id[i_pt];
                if (!table->hasKey(pt_id)) continue;
                active_pt_list.push_back(pt_id);
                active_eqs_id_list.push_back(active_entry_counter++);
            }
            local_table.addVariableDoFs(mesh_id, active_pt_list);
            IEquationIdStorage *local_storage = local_table.getPointEquationIdTable(i_var, mesh_id);
            for (std::size_t j=0; j<active_pt_list.size(); j++) {
                local_storage->set(active_pt_list[j], active_eqs_id_list[j]);
            }
        }
    } else { // by point
        std::size_t active_entry_counter = 0;
        std::vector<std::vector<std::size_t> > var_active_pt_list(list_var.size());
        std::vector<std::vector<std::size_t> > var_active_eqs_id_list(list_var.size());
        for (std::size_t i_pt=0; i_pt<n_pt; i_pt++) {
            const std::size_t pt_id = list_pt_id[i_pt];
            for (std::size_t i_var=0; i_var<list_var.size(); i_var++) {
                const std::size_t var_id = list_var[i_var];
                const IEquationIdStorage* table = getPointEquationIdTable(var_id, mesh_id);
                if (!table->hasKey(pt_id)) continue;
                var_active_pt_list[i_var].push_back(pt_id);
                var_active_eqs_id_list[i_var].push_back(active_entry_counter++);
            }
        }

        for (std::size_t i_var=0; i_var<var_active_pt_list.size(); i_var++) {
            local_table.addVariableDoFs(mesh_id, var_active_pt_list[i_var]);
            IEquationIdStorage *local_storage = local_table.getPointEquationIdTable(i_var, mesh_id);
            for (std::size_t j=0; j<var_active_pt_list[i_var].size(); j++) {
                local_storage->set(var_active_pt_list[i_var][j], var_active_eqs_id_list[i_var][j]);
            }
        }
    }

    local_table.setNumberingType(_local_numbering_type);
    //local_table.construct();
}

void DofEquationIdTable::mapDoF(std::size_t eqs_id, std::size_t &var_id, std::size_t &mesh_id, std::size_t &pt_id) const
{
    var_id = IEquationIdStorage::index_npos;
    mesh_id = IEquationIdStorage::index_npos;
    pt_id = IEquationIdStorage::index_npos;
    for (std::size_t i=0; i<_map_var2dof.size(); i++) {
        const std::map<std::size_t, IEquationIdStorage*> &obj = _map_var2dof[i];
        for (std::map<std::size_t, IEquationIdStorage*>::const_iterator itr=obj.begin(); itr!=obj.end(); ++itr) {
            const IEquationIdStorage* pt2dof = itr->second;
            if (!pt2dof->hasValue(eqs_id)) continue;
            pt_id = pt2dof->key(eqs_id);
            var_id = i;
            mesh_id = itr->first;
            return;
        }
    }
}

const IEquationIdStorage* DofEquationIdTable::getPointEquationIdTable(std::size_t var_id, std::size_t mesh_id) const
{
    std::map<std::size_t, IEquationIdStorage*>::const_iterator itr = _map_var2dof[var_id].find(mesh_id);
    if (itr!=_map_var2dof[var_id].end()) {
        return itr->second;
    } else {
        return 0;
    }
}

void DofEquationIdTable::construct(long offset)
{
    if (_numbering_type==DofNumberingType::BY_VARIABLE) {
        //order by dof
        long eqs_id = offset;
        for (std::size_t i=0; i<_map_var2dof.size(); i++) {
            std::map<std::size_t, IEquationIdStorage*> &obj = _map_var2dof[i];
            for (std::map<std::size_t, IEquationIdStorage*>::iterator itr = obj.begin(); itr!=obj.end(); ++itr) {
                IEquationIdStorage* pt2eq = itr->second;
                eqs_id = pt2eq->setAll(eqs_id);
            }
        }
        _total_dofs_without_ghosts = eqs_id; //TODO
        _total_dofs = eqs_id;
    } else {
        //order by discrete points
        long eqs_id = offset;
        if (_is_seq_address) {
            assert(_ghost_pt.size()==0);
            for (std::map<std::size_t, std::vector<std::size_t> >::iterator itr=_map_msh2var.begin(); itr!=_map_msh2var.end(); ++itr) {
                std::size_t mesh_id = itr->first;
                std::vector<std::size_t> &list_var = itr->second;
                for (std::size_t i=0; i<list_var.size(); i++) {
                    IEquationIdStorage* pt2eq = getPointEquationIdTable(list_var[i], mesh_id);
                    eqs_id = offset + i;
                    eqs_id = pt2eq->setAll(eqs_id, list_var.size());
                }
            }
            _total_dofs_without_ghosts = eqs_id; //TODO
            _total_dofs = eqs_id;
        } else { //random address
            for (std::map<std::size_t, std::vector<std::size_t> >::iterator itr=_map_msh2var.begin(); itr!=_map_msh2var.end(); ++itr) {
                std::size_t mesh_id = itr->first;
                std::vector<std::size_t> &list_var = itr->second;
                // get range of point id
                std::size_t i_min = std::numeric_limits<std::size_t>::max();
                std::size_t i_max = 0;
                for (std::size_t i=0; i<list_var.size(); i++) {
                    IEquationIdStorage* pt2eq = getPointEquationIdTable(list_var[i], mesh_id);
                    std::size_t tmp_min, tmp_max;
                    pt2eq->key_range(tmp_min, tmp_max);
                    i_min = std::min(i_min, tmp_min);
                    i_max = std::max(i_max, tmp_max);
                }
                // first, real nodes
                for (std::size_t i=i_min; i<i_max+1; i++) {
                    if (isGhostPoint(mesh_id, i)) continue; //skip ghost
                    for (std::size_t j=0; j<list_var.size(); j++) {
                        IEquationIdStorage* pt2eq = getPointEquationIdTable(list_var[j], mesh_id);
                        if (pt2eq->hasKey(i) && pt2eq->isActive(i)) {
                            pt2eq->set(i, eqs_id++);
                        }
                    }
                }
                _total_dofs_without_ghosts = eqs_id;
                // next ghost nodes
                for (std::size_t i=i_min; i<i_max+1; i++) {
                    if (!isGhostPoint(mesh_id, i)) continue; //skip real
                    for (std::size_t j=0; j<list_var.size(); j++) {
                        IEquationIdStorage* pt2eq = getPointEquationIdTable(list_var[j], mesh_id);
                        if (pt2eq->hasKey(i) && pt2eq->isActive(i)) {
                            pt2eq->set(i, eqs_id++);
                        }
                    }
                }
            }
            _total_dofs = eqs_id;
        }
    }
}



} //end
