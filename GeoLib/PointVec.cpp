/**
 * \file
 * \author Thomas Fischeror
 * \date   2010-06-11
 * \brief  Implementation of the PointVec class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <numeric>

#include "BaseLib/Logging.h"

#include "PointVec.h"

#include "MathLib/MathTools.h"

namespace GeoLib
{
PointVec::PointVec(
    const std::string& name, std::unique_ptr<std::vector<Point*>> points,
    std::unique_ptr<std::map<std::string, std::size_t>> name_id_map,
    PointType type, double rel_eps)
    : TemplateVec<Point>(name, std::move(points), std::move(name_id_map)),
      type_(type),
      aabb_(data_vec_->begin(), data_vec_->end()),
      rel_eps_(rel_eps * std::sqrt(MathLib::sqrDist(aabb_.getMinPoint(),
                                                    aabb_.getMaxPoint()))),
      oct_tree_(OctTree<GeoLib::Point, 16>::createOctTree(
          aabb_.getMinPoint(), aabb_.getMaxPoint(), rel_eps_))
{
    assert(data_vec_);
    std::size_t const number_of_all_input_pnts(data_vec_->size());

    // correct the ids if necessary
    for (std::size_t k(0); k < data_vec_->size(); ++k)
    {
        if ((*data_vec_)[k]->getID() == std::numeric_limits<std::size_t>::max())
        {
            (*data_vec_)[k]->setID(k);
        }
    }

    std::vector<std::size_t> rm_pos;
    // add all points in the oct tree in order to make them unique
    pnt_id_map_.resize(number_of_all_input_pnts);
    std::iota(pnt_id_map_.begin(), pnt_id_map_.end(), 0);
    GeoLib::Point* ret_pnt(nullptr);
    for (std::size_t k(0); k < data_vec_->size(); ++k)
    {
        GeoLib::Point* const pnt((*data_vec_)[k]);
        if (!oct_tree_->addPoint(pnt, ret_pnt))
        {
            assert(ret_pnt != nullptr);
            pnt_id_map_[pnt->getID()] = ret_pnt->getID();
            rm_pos.push_back(k);
            delete (*data_vec_)[k];
            (*data_vec_)[k] = nullptr;
        }
        else
        {
            pnt_id_map_[k] = pnt->getID();
        }
    }

    auto const data_vec_end =
        std::remove(data_vec_->begin(), data_vec_->end(), nullptr);
    data_vec_->erase(data_vec_end, data_vec_->end());

    // decrement the ids according to the number of removed points (==k) before
    // the j-th point (positions of removed points are stored in the vector
    // rm_pos)
    for (std::size_t k(1); k < rm_pos.size(); ++k)
    {
        // decrement the ids in the interval [rm_pos[k-1]+1, rm_pos[k])
        for (std::size_t j(rm_pos[k - 1] + 1); j < rm_pos[k]; ++j)
        {
            pnt_id_map_[j] -= k;
        }
    }
    // decrement the ids from rm_pos.back()+1 until the end of pnt_id_map_
    if (!rm_pos.empty())
    {
        for (std::size_t j(rm_pos.back() + 1); j < pnt_id_map_.size(); ++j)
        {
            pnt_id_map_[j] -= rm_pos.size();
        }
    }
    // decrement the ids within the pnt_id_map_ at positions of the removed
    // points
    for (std::size_t k(1); k < rm_pos.size(); ++k)
    {
        std::size_t cnt(0);
        for (cnt = 0;
             cnt < rm_pos.size() && pnt_id_map_[rm_pos[k]] > rm_pos[cnt];)
        {
            cnt++;
        }
        pnt_id_map_[rm_pos[k]] -= cnt;
    }

    // set value of the point id to the position of the point within data_vec_
    for (std::size_t k(0); k < data_vec_->size(); ++k)
    {
        (*data_vec_)[k]->setID(k);
    }

    if (number_of_all_input_pnts > data_vec_->size())
    {
        WARN("PointVec::PointVec(): there are {:d} double points.",
             number_of_all_input_pnts - data_vec_->size());
    }

    correctNameIDMapping();
    // create the inverse mapping
    id_to_name_map_.resize(data_vec_->size());
    // fetch the names from the name id map
    for (auto p : *name_id_map_)
    {
        if (p.second >= id_to_name_map_.size())
        {
            continue;
        }
        id_to_name_map_[p.second] = p.first;
    }
}

std::size_t PointVec::push_back(Point* pnt)
{
    pnt_id_map_.push_back(uniqueInsert(pnt));
    id_to_name_map_.emplace_back("");
    return pnt_id_map_[pnt_id_map_.size() - 1];
}

void PointVec::push_back(Point* pnt, std::string const* const name)
{
    if (name == nullptr)
    {
        pnt_id_map_.push_back(uniqueInsert(pnt));
        id_to_name_map_.emplace_back("");
        return;
    }

    std::map<std::string, std::size_t>::const_iterator it(
        name_id_map_->find(*name));
    if (it != name_id_map_->end())
    {
        id_to_name_map_.emplace_back("");
        WARN("PointVec::push_back(): two points share the name {:s}.",
             name->c_str());
        return;
    }

    std::size_t id(uniqueInsert(pnt));
    pnt_id_map_.push_back(id);
    (*name_id_map_)[*name] = id;
    id_to_name_map_.push_back(*name);
}

std::size_t PointVec::uniqueInsert(Point* pnt)
{
    GeoLib::Point* ret_pnt(nullptr);
    if (oct_tree_->addPoint(pnt, ret_pnt))
    {
        // set value of the point id to the position of the point within
        // data_vec_
        pnt->setID(data_vec_->size());
        data_vec_->push_back(pnt);
        return data_vec_->size() - 1;
    }

    // pnt is outside of OctTree object
    if (ret_pnt == nullptr)
    {
        // update the axis aligned bounding box
        aabb_.update(*pnt);
        // recreate the (enlarged) OctTree
        oct_tree_.reset(GeoLib::OctTree<GeoLib::Point, 16>::createOctTree(
            aabb_.getMinPoint(), aabb_.getMaxPoint(), rel_eps_));
        // add all points that are already in the data_vec_
        for (std::size_t k(0); k < data_vec_->size(); ++k)
        {
            GeoLib::Point* const p((*data_vec_)[k]);
            oct_tree_->addPoint(p, ret_pnt);
        }
        // add the new point
        ret_pnt = nullptr;
        oct_tree_->addPoint(pnt, ret_pnt);
        // set value of the point id to the position of the point within
        // data_vec_
        pnt->setID(data_vec_->size());
        data_vec_->push_back(pnt);
        return data_vec_->size() - 1;
    }

    delete pnt;
    return ret_pnt->getID();
}

void PointVec::correctNameIDMapping()
{
    // create mapping id -> name using the std::vector id_names
    std::vector<std::string> id_names(pnt_id_map_.size(), std::string(""));
    for (auto& id_name_pair : *name_id_map_)
    {
        id_names[id_name_pair.second] = id_name_pair.first;
    }

    for (auto it = name_id_map_->begin(); it != name_id_map_->end();)
    {
        // extract the id associated with the name
        const std::size_t id(it->second);

        if (pnt_id_map_[id] == id)
        {
            ++it;
            continue;
        }

        if (pnt_id_map_[pnt_id_map_[id]] == pnt_id_map_[id])
        {
            if (id_names[pnt_id_map_[id]].length() != 0)
            {
                // point has already a name, erase the second occurrence
                it = name_id_map_->erase(it);
            }
            else
            {
                // until now the point has not a name assign the second
                // occurrence the correct id
                it->second = pnt_id_map_[id];
                ++it;
            }
        }
        else
        {
            it->second = pnt_id_map_[id];  // update id associated to the name
            ++it;
        }
    }
}

std::string const& PointVec::getItemNameByID(std::size_t id) const
{
    return id_to_name_map_[id];
}

void PointVec::setNameForElement(std::size_t id, std::string const& name)
{
    TemplateVec::setNameForElement(id, name);
    id_to_name_map_[id] = name;
}

void PointVec::resetInternalDataStructures()
{
    MathLib::Point3d const& min(aabb_.getMinPoint());
    MathLib::Point3d const& max(aabb_.getMaxPoint());
    double const rel_eps(rel_eps_ / std::sqrt(MathLib::sqrDist(min, max)));

    aabb_ = GeoLib::AABB(data_vec_->begin(), data_vec_->end());

    rel_eps_ = rel_eps * std::sqrt(MathLib::sqrDist(aabb_.getMinPoint(),
                                                    aabb_.getMaxPoint()));

    oct_tree_.reset(OctTree<GeoLib::Point, 16>::createOctTree(
        aabb_.getMinPoint(), aabb_.getMaxPoint(), rel_eps_));

    GeoLib::Point* ret_pnt(nullptr);
    for (auto const p : *data_vec_)
    {
        oct_tree_->addPoint(p, ret_pnt);
    }
}

}  // namespace GeoLib
