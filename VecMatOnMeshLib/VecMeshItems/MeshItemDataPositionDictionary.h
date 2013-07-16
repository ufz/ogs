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


#ifndef MESHITEMDATAPOSITIONDICTIONARY_H_
#define MESHITEMDATAPOSITIONDICTIONARY_H_

#include <iostream>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>

namespace VecMatOnMeshLib
{

struct MeshitemDataPosition
{
    // Three ids for a location in a mesh
    std::size_t mesh_id;
    MeshItemType::type item_type;
    std::size_t mesh_item_id;

    // Physical component
    std::size_t comp_id;

    // Position in global matrix or vector
    std::size_t global_index;

    MeshitemDataPosition(std::size_t mesh_id, MeshItemType::type item_type, std::size_t mesh_item_id, std::size_t comp_id, std::size_t global_index)
    : mesh_id(mesh_id), item_type(item_type), mesh_item_id(mesh_item_id), comp_id(comp_id), global_index(global_index)
    {}

    void print() const
    {
        std::cout << mesh_id << ", " << item_type << ", " << mesh_item_id << ", " << comp_id << ", " << global_index << "\n";
    }
};

struct mesh_item_ID {};
struct mesh_item_comp_ID {};
struct comp_ID {};
struct ByGlobalIndex {};

struct mesh_item_key : public boost::multi_index::composite_key<
    MeshitemDataPosition,
    BOOST_MULTI_INDEX_MEMBER(MeshitemDataPosition,std::size_t,mesh_id),
    BOOST_MULTI_INDEX_MEMBER(MeshitemDataPosition,MeshItemType::type,item_type),
    BOOST_MULTI_INDEX_MEMBER(MeshitemDataPosition,std::size_t,mesh_item_id)
    >{};

struct mesh_item_comp_key: public boost::multi_index::composite_key<
    MeshitemDataPosition,
    BOOST_MULTI_INDEX_MEMBER(MeshitemDataPosition,std::size_t,mesh_id),
    BOOST_MULTI_INDEX_MEMBER(MeshitemDataPosition,MeshItemType::type,item_type),
    BOOST_MULTI_INDEX_MEMBER(MeshitemDataPosition,std::size_t,mesh_item_id),
    BOOST_MULTI_INDEX_MEMBER(MeshitemDataPosition,std::size_t,comp_id)
    >{};

typedef boost::multi_index::multi_index_container<
        MeshitemDataPosition,
        boost::multi_index::indexed_by
        <
            boost::multi_index::ordered_unique
            <
                boost::multi_index::tag<mesh_item_comp_ID>, mesh_item_comp_key
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<mesh_item_ID>, mesh_item_key
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<comp_ID>,
                boost::multi_index::member<MeshitemDataPosition, std::size_t, &MeshitemDataPosition::comp_id>
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByGlobalIndex>,
                boost::multi_index::member<MeshitemDataPosition, std::size_t, &MeshitemDataPosition::global_index>
            >
        >
    > MeshitemDataPositionDictionary;

} // VecMatOnMeshLib

#endif /* MESHITEMDATAPOSITIONDICTIONARY_H_ */
