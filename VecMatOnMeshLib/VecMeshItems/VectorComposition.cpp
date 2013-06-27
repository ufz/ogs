/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include "VectorComposition.h"

#include <numeric>
#include <iostream>
#include <boost/range/algorithm/for_each.hpp>

namespace VecMatOnMeshLib
{

VectorComposition::VectorComposition(const std::vector<ComponentDistribution*> &vec_comp_dis, OrderingType::type numbering)
: _vec_comp_dis(vec_comp_dis), _ordering_type(numbering)
{
    _n_data = std::accumulate(vec_comp_dis.begin(), vec_comp_dis.end(),
                                        0u,
                                        [](std::size_t sum, const ComponentDistribution* comp_dis)
                                        {
                                            return sum+comp_dis->getNMeshItems();
                                        }
                                       );

    // construct dict (and here we number data_id by component type)
    std::size_t data_id = 0;
    for (auto itrComp = _vec_comp_dis.begin(); itrComp!=_vec_comp_dis.end(); ++itrComp) {
        auto comp_id = distance(_vec_comp_dis.begin(),itrComp);
        for (unsigned i=0; i<(*itrComp)->getNMeshes(); i++) {
            auto mesh_items = (*itrComp)->getMeshItems(i);
            std::size_t mesh_id = mesh_items.getMeshID();
            _vec_meshIDs.insert(mesh_id);
            // mesh items are ordered first by node, cell, ....
            for (std::size_t j=0; j<mesh_items.getNNodes(); j++) {
                _dict.insert(MeshitemDataPosition(mesh_id, MeshItemType::Node, j, comp_id, data_id++));
            }
            for (std::size_t j=0; j<mesh_items.getNElements(); j++) {
                _dict.insert(MeshitemDataPosition(mesh_id, MeshItemType::Cell, j, comp_id, data_id++));
            }
        }
    }

    if (numbering!=OrderingType::BY_COMPONENT_TYPE)
        numberingByMeshItems();
}


void VectorComposition::numberingByComponent(std::size_t offset)
{
    std::size_t data_id = offset;

    auto &m = _dict.get<comp_ID>(); // view as sorted by comp ID
    for (auto itr_mesh_item=m.begin(); itr_mesh_item!=m.end(); ++itr_mesh_item) {
        MeshitemDataPosition pos = *itr_mesh_item;
        pos.data_id = data_id++;
        m.replace(itr_mesh_item, pos);
    }
//    std::size_t data_id = offset;
//    for (auto itrComp = _vec_comp_dis.begin(); itrComp!=_vec_comp_dis.end(); ++itrComp) {
//        for (unsigned i=0; i<(*itrComp)->getNMeshes(); i++) {
//            auto mesh_items = (*itrComp)->getMeshItems(i);
//            for (std::size_t j=0; j<mesh_items.getNTotalItems(); j++) {
//                _dict.insert(MeshitemDataPosition(mesh_items.getMeshID(), j, (*itrComp)->getComponentID(), data_id));
//                data_id++;
//            }
//        }
//    }
}

void VectorComposition::numberingByMeshItems(std::size_t offset)
{
    std::size_t data_id = offset;

    auto &m = _dict.get<mesh_item_ID>(); // view as sorted by mesh item
    for (auto itr_mesh_item=m.begin(); itr_mesh_item!=m.end(); ++itr_mesh_item) {
        MeshitemDataPosition pos = *itr_mesh_item;
        pos.data_id = data_id++;
        m.replace(itr_mesh_item, pos);
    }
}


std::size_t VectorComposition::getDataID(const MeshItem &pos, unsigned compID) const
{
    auto &m = _dict.get<mesh_item_comp_ID>();
    auto itr = m.find(boost::make_tuple(pos.mesh_id, pos.item_type, pos.item_id, compID));
    return itr!=m.end() ? itr->data_id : -1;
}

std::vector<std::size_t> VectorComposition::getComponentIDs(const MeshItem &pos) const
{
    auto &m = _dict.get<mesh_item_ID>();
    auto p = m.equal_range(boost::make_tuple(pos.mesh_id, pos.item_type, pos.item_id));
    std::vector<std::size_t> vec_compID;
    for (auto itr=p.first; itr!=p.second; ++itr)
        vec_compID.push_back(itr->comp_id);
    return vec_compID;
}

std::vector<std::size_t> VectorComposition::getDataIDList(const MeshItem &pos) const
{
    auto &m = _dict.get<mesh_item_ID>();
    auto p = m.equal_range(boost::make_tuple(pos.mesh_id, pos.item_type, pos.item_id));
    std::vector<std::size_t> vec_dataID;
    for (auto itr=p.first; itr!=p.second; ++itr)
        vec_dataID.push_back(itr->data_id);
    return vec_dataID;
}

std::vector<std::size_t> VectorComposition::getDataIDList(const std::vector<MeshItem> &vec_pos, OrderingType::type list_numbering) const
{
    MeshitemDataPositionDictionary sub_dict;
    {
        auto &m = _dict.get<mesh_item_ID>();
        for (auto &item : vec_pos) {
            auto p = m.equal_range(boost::make_tuple(item.mesh_id, item.item_type, item.item_id));
            for (auto itr=p.first; itr!=p.second; ++itr)
                sub_dict.insert(*itr);
        }
    }

    std::vector<std::size_t> vec_dataID;
    if (list_numbering==OrderingType::BY_MESH_ITEM_ID) {
        auto &m = sub_dict.get<mesh_item_ID>();
        for (auto itr_mesh_item=m.begin(); itr_mesh_item!=m.end(); ++itr_mesh_item) {
            vec_dataID.push_back(itr_mesh_item->data_id);
        }
    } else {
        auto &m = sub_dict.get<comp_ID>();
        for (auto itr_mesh_item=m.begin(); itr_mesh_item!=m.end(); ++itr_mesh_item) {
            vec_dataID.push_back(itr_mesh_item->data_id);
        }
    }
    return vec_dataID;
}

void disp(const MeshitemDataPosition &dat)
{
    dat.print();
}

void VectorComposition::print()
{

    boost::for_each(_dict, disp);
}

} // VecMatOnMeshLib

