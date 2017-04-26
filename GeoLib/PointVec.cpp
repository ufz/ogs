/**
 * \file
 * \author Thomas Fischeror
 * \date   2010-06-11
 * \brief  Implementation of the PointVec class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <numeric>

#include <logog/include/logog.hpp>

#include "PointVec.h"

#include "MathLib/MathTools.h"

namespace GeoLib
{
PointVec::PointVec(const std::string& name,
                   std::unique_ptr<std::vector<Point*>> points,
                   std::map<std::string, std::size_t>* name_id_map,
                   PointType type, double rel_eps)
    : TemplateVec<Point>(name, std::move(points), name_id_map),
      _type(type),
      _aabb(_data_vec->begin(), _data_vec->end()),
      _rel_eps(rel_eps * std::sqrt(MathLib::sqrDist(_aabb.getMinPoint(),
                                                    _aabb.getMaxPoint()))),
      _oct_tree(OctTree<GeoLib::Point, 16>::createOctTree(
          _aabb.getMinPoint(), _aabb.getMaxPoint(), _rel_eps))
{
    assert(_data_vec);
    std::size_t const number_of_all_input_pnts(_data_vec->size());

    // correct the ids if necessary
    for (std::size_t k(0); k < _data_vec->size(); ++k)
    {
        if ((*_data_vec)[k]->getID() == std::numeric_limits<std::size_t>::max())
        {
            (*_data_vec)[k]->setID(k);
        }
    }

    std::vector<std::size_t> rm_pos;
    // add all points in the oct tree in order to make them unique
    _pnt_id_map.resize(number_of_all_input_pnts);
    std::iota(_pnt_id_map.begin(), _pnt_id_map.end(), 0);
    GeoLib::Point* ret_pnt(nullptr);
    for (std::size_t k(0); k < _data_vec->size(); ++k)
    {
        GeoLib::Point* const pnt((*_data_vec)[k]);
        if (!_oct_tree->addPoint(pnt, ret_pnt))
        {
            assert(ret_pnt != nullptr);
            _pnt_id_map[pnt->getID()] = ret_pnt->getID();
            rm_pos.push_back(k);
            delete (*_data_vec)[k];
            (*_data_vec)[k] = nullptr;
        }
        else
        {
            _pnt_id_map[k] = pnt->getID();
        }
    }

    auto const data_vec_end =
        std::remove(_data_vec->begin(), _data_vec->end(), nullptr);
    _data_vec->erase(data_vec_end, _data_vec->end());

    // decrement the ids according to the number of removed points (==k) before
    // the j-th point (positions of removed points are stored in the vector
    // rm_pos)
    for (std::size_t k(1); k < rm_pos.size(); ++k)
    {
        // decrement the ids in the interval [rm_pos[k-1]+1, rm_pos[k])
        for (std::size_t j(rm_pos[k - 1] + 1); j < rm_pos[k]; ++j)
        {
            _pnt_id_map[j] -= k;
        }
    }
    // decrement the ids from rm_pos.back()+1 until the end of _pnt_id_map
    if (!rm_pos.empty())
    {
        for (std::size_t j(rm_pos.back() + 1); j < _pnt_id_map.size(); ++j)
        {
            _pnt_id_map[j] -= rm_pos.size();
        }
    }
    // decrement the ids within the _pnt_id_map at positions of the removed
    // points
    for (std::size_t k(1); k < rm_pos.size(); ++k)
    {
        std::size_t cnt(0);
        for (cnt = 0;
             cnt < rm_pos.size() && _pnt_id_map[rm_pos[k]] > rm_pos[cnt];)
            cnt++;
        _pnt_id_map[rm_pos[k]] -= cnt;
    }

    // set value of the point id to the position of the point within _data_vec
    for (std::size_t k(0); k < _data_vec->size(); ++k)
        (*_data_vec)[k]->setID(k);

    if (number_of_all_input_pnts > _data_vec->size())
        WARN("PointVec::PointVec(): there are %d double points.",
             number_of_all_input_pnts - _data_vec->size());

    correctNameIDMapping();
    // create the inverse mapping
    _id_to_name_map.resize(_data_vec->size());
    // fetch the names from the name id map
    for (auto p : *_name_id_map)
    {
        if (p.second >= _id_to_name_map.size()) continue;
        _id_to_name_map[p.second] = p.first;
    }
}

std::size_t PointVec::push_back(Point* pnt)
{
    _pnt_id_map.push_back(uniqueInsert(pnt));
    _id_to_name_map.emplace_back("");
    return _pnt_id_map[_pnt_id_map.size() - 1];
}

void PointVec::push_back(Point* pnt, std::string const* const name)
{
    if (name == nullptr)
    {
        _pnt_id_map.push_back(uniqueInsert(pnt));
        _id_to_name_map.emplace_back("");
        return;
    }

    std::map<std::string, std::size_t>::const_iterator it(
        _name_id_map->find(*name));
    if (it != _name_id_map->end())
    {
        _id_to_name_map.emplace_back("");
        WARN("PointVec::push_back(): two points share the name %s.",
             name->c_str());
        return;
    }

    std::size_t id(uniqueInsert(pnt));
    _pnt_id_map.push_back(id);
    (*_name_id_map)[*name] = id;
    _id_to_name_map.push_back(*name);
}

std::size_t PointVec::uniqueInsert(Point* pnt)
{
    GeoLib::Point* ret_pnt(nullptr);
    if (_oct_tree->addPoint(pnt, ret_pnt))
    {
        // set value of the point id to the position of the point within
        // _data_vec
        pnt->setID(_data_vec->size());
        _data_vec->push_back(pnt);
        return _data_vec->size() - 1;
    }

    // pnt is outside of OctTree object
    if (ret_pnt == nullptr)
    {
        // update the axis aligned bounding box
        _aabb.update(*pnt);
        // recreate the (enlarged) OctTree
        _oct_tree.reset(GeoLib::OctTree<GeoLib::Point, 16>::createOctTree(
            _aabb.getMinPoint(), _aabb.getMaxPoint(), _rel_eps));
        // add all points that are already in the _data_vec
        for (std::size_t k(0); k < _data_vec->size(); ++k)
        {
            GeoLib::Point* const p((*_data_vec)[k]);
            _oct_tree->addPoint(p, ret_pnt);
        }
        // add the new point
        ret_pnt = nullptr;
        _oct_tree->addPoint(pnt, ret_pnt);
        // set value of the point id to the position of the point within
        // _data_vec
        pnt->setID(_data_vec->size());
        _data_vec->push_back(pnt);
        return _data_vec->size() - 1;
    }
    else
    {
        delete pnt;
        return ret_pnt->getID();
    }
}

void PointVec::correctNameIDMapping()
{
    // create mapping id -> name using the std::vector id_names
    std::vector<std::string> id_names(_pnt_id_map.size(), std::string(""));
    for (auto& id_name_pair : *_name_id_map)
    {
        id_names[id_name_pair.second] = id_name_pair.first;
    }

    for (auto it = _name_id_map->begin(); it != _name_id_map->end();)
    {
        // extract the id associated with the name
        const std::size_t id(it->second);

        if (_pnt_id_map[id] == id)
        {
            ++it;
            continue;
        }

        if (_pnt_id_map[_pnt_id_map[id]] == _pnt_id_map[id])
        {
            if (id_names[_pnt_id_map[id]].length() != 0)
            {
                // point has already a name, erase the second occurrence
                it = _name_id_map->erase(it);
            }
            else
            {
                // until now the point has not a name assign the second
                // occurrence the correct id
                it->second = _pnt_id_map[id];
                ++it;
            }
        }
        else
        {
            it->second = _pnt_id_map[id];  // update id associated to the name
            ++it;
        }
    }
}

std::string const& PointVec::getItemNameByID(std::size_t id) const
{
    return _id_to_name_map[id];
}

void PointVec::setNameForElement(std::size_t id, std::string const& name)
{
    TemplateVec::setNameForElement(id, name);
    _id_to_name_map[id] = name;
}

void PointVec::resetInternalDataStructures()
{
    MathLib::Point3d const& min(_aabb.getMinPoint());
    MathLib::Point3d const& max(_aabb.getMaxPoint());
    double const rel_eps(_rel_eps / std::sqrt(MathLib::sqrDist(min, max)));

    _aabb = GeoLib::AABB(_data_vec->begin(), _data_vec->end());

    _rel_eps = rel_eps * std::sqrt(MathLib::sqrDist(_aabb.getMinPoint(),
                                                    _aabb.getMaxPoint()));

    _oct_tree.reset(OctTree<GeoLib::Point, 16>::createOctTree(
        _aabb.getMinPoint(), _aabb.getMaxPoint(), _rel_eps));

    GeoLib::Point* ret_pnt(nullptr);
    for (auto const p : *_data_vec)
    {
        _oct_tree->addPoint(p, ret_pnt);
    }
}

}  // end namespace
