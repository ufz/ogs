/**
 * \date   2015-06-12
 * \brief  Implementation of the OctTree class.
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

namespace GeoLib
{
template <typename POINT, std::size_t MAX_POINTS>
OctTree<POINT, MAX_POINTS>* OctTree<POINT, MAX_POINTS>::createOctTree(
    Eigen::Vector3d ll, Eigen::Vector3d ur, double eps)
{
    // compute an axis aligned cube around the points ll and ur
    const double dx(ur[0] - ll[0]);
    const double dy(ur[1] - ll[1]);
    const double dz(ur[2] - ll[2]);
    if (dx >= dy && dx >= dz)
    {
        ll[1] -= (dx - dy) / 2.0;
        ur[1] += (dx - dy) / 2.0;
        ll[2] -= (dx - dz) / 2.0;
        ur[2] += (dx - dz) / 2.0;
    }
    else
    {
        if (dy >= dx && dy >= dz)
        {
            ll[0] -= (dy - dx) / 2.0;
            ur[0] += (dy - dx) / 2.0;
            ll[2] -= (dy - dz) / 2.0;
            ur[2] += (dy - dz) / 2.0;
        }
        else
        {
            ll[0] -= (dz - dx) / 2.0;
            ur[0] += (dz - dx) / 2.0;
            ll[1] -= (dz - dy) / 2.0;
            ur[1] += (dz - dy) / 2.0;
        }
    }
    if (eps == 0.0)
    {
        eps = std::numeric_limits<double>::epsilon();
    }
    for (std::size_t k(0); k < 3; ++k)
    {
        if (ur[k] - ll[k] > 0.0)
        {
            ur[k] += (ur[k] - ll[k]) * 1e-6;
        }
        else
        {
            ur[k] += eps;
        }
    }
    Eigen::Vector3d lower_left{ll[0], ll[1], ll[2]};
    Eigen::Vector3d upper_right{ur[0], ur[1], ur[2]};
    return new OctTree<POINT, MAX_POINTS>(lower_left, upper_right, eps);
}

template <typename POINT, std::size_t MAX_POINTS>
OctTree<POINT, MAX_POINTS>::~OctTree()
{
    for (auto c : _children)
    {
        delete c;
    }
}

template <typename POINT, std::size_t MAX_POINTS>
bool OctTree<POINT, MAX_POINTS>::addPoint(POINT* pnt, POINT*& ret_pnt)
{
    // first do a range query using a epsilon box around the point pnt
    std::vector<POINT*> query_pnts;
    Eigen::Vector3d const min =
        Eigen::Map<Eigen::Vector3d>(pnt->getCoords()).array() - _eps;
    Eigen::Vector3d const max =
        Eigen::Map<Eigen::Vector3d>(pnt->getCoords()).array() + _eps;
    getPointsInRange(min, max, query_pnts);
    auto const it = std::find_if(
        query_pnts.begin(), query_pnts.end(),
        [pnt, this](auto const* p)
        {
            return (Eigen::Map<Eigen::Vector3d const>(p->getCoords()) -
                    Eigen::Map<Eigen::Vector3d const>(pnt->getCoords()))
                       .squaredNorm() < _eps * _eps;
        });
    if (it != query_pnts.end())
    {
        ret_pnt = *it;
        return false;
    }

    // the point pnt is not yet in the OctTree
    if (isOutside(pnt))
    {
        ret_pnt = nullptr;
        return false;
    }

    // at this place it holds true that the point is within [_ll, _ur]
    if (!_is_leaf)
    {
        for (auto c : _children)
        {
            if (c->addPoint_(pnt, ret_pnt))
            {
                return true;
            }
            if (ret_pnt != nullptr)
            {
                return false;
            }
        }
    }

    ret_pnt = pnt;

    if (_pnts.size() < MAX_POINTS)
    {
        _pnts.push_back(pnt);
    }
    else
    {  // i.e. _pnts.size () == MAX_POINTS
        splitNode(pnt);
        _pnts.clear();
    }
    return true;
}

template <typename POINT, std::size_t MAX_POINTS>
template <typename T>
void OctTree<POINT, MAX_POINTS>::getPointsInRange(
    T const& min, T const& max, std::vector<POINT*>& pnts) const
{
    if (_ur[0] < min[0] || _ur[1] < min[1] || _ur[2] < min[2])
    {
        return;
    }

    if (max[0] < _ll[0] || max[1] < _ll[1] || max[2] < _ll[2])
    {
        return;
    }

    if (_is_leaf)
    {
        std::copy_if(_pnts.begin(), _pnts.end(), std::back_inserter(pnts),
                     [&min, &max](auto const* p)
                     {
                         return (min[0] <= (*p)[0] && (*p)[0] < max[0] &&
                                 min[1] <= (*p)[1] && (*p)[1] < max[1] &&
                                 min[2] <= (*p)[2] && (*p)[2] < max[2]);
                     });
    }
    else
    {
        for (std::size_t k(0); k < 8; k++)
        {
            _children[k]->getPointsInRange(min, max, pnts);
        }
    }
}

template <typename POINT, std::size_t MAX_POINTS>
OctTree<POINT, MAX_POINTS>::OctTree(Eigen::Vector3d const& ll,
                                    Eigen::Vector3d const& ur, double eps)
    : _ll(ll), _ur(ur), _is_leaf(true), _eps(eps)
{
    _children.fill(nullptr);
}

template <typename POINT, std::size_t MAX_POINTS>
bool OctTree<POINT, MAX_POINTS>::addPoint_(POINT* pnt, POINT*& ret_pnt)
{
    if (isOutside(pnt))
    {
        ret_pnt = nullptr;
        return false;
    }

    // at this place it holds true that the point is within [_ll, _ur]
    if (!_is_leaf)
    {
        for (auto c : _children)
        {
            if (c->addPoint_(pnt, ret_pnt))
            {
                return true;
            }
            if (ret_pnt != nullptr)
            {
                return false;
            }
        }
    }

    ret_pnt = pnt;
    if (_pnts.size() < MAX_POINTS)
    {
        _pnts.push_back(pnt);
    }
    else
    {  // i.e. _pnts.size () == MAX_POINTS
        splitNode(pnt);
        _pnts.clear();
    }
    return true;
}

template <typename POINT, std::size_t MAX_POINTS>
bool OctTree<POINT, MAX_POINTS>::addPointToChild(POINT* pnt)
{
    if (isOutside(pnt))
    {
        return false;
    }

    if (_pnts.size() < MAX_POINTS)
    {
        _pnts.push_back(pnt);
    }
    else
    {  // i.e. _pnts.size () == MAX_POINTS
        splitNode(pnt);
        _pnts.clear();
    }
    return true;
}

template <typename POINT, std::size_t MAX_POINTS>
void OctTree<POINT, MAX_POINTS>::splitNode(POINT* pnt)
{
    const double x_mid((_ur[0] + _ll[0]) / 2.0);
    const double y_mid((_ur[1] + _ll[1]) / 2.0);
    const double z_mid((_ur[2] + _ll[2]) / 2.0);
    Eigen::Vector3d p0{x_mid, y_mid, _ll[2]};
    Eigen::Vector3d p1{_ur[0], _ur[1], z_mid};

    // create child NEL
    _children[static_cast<std::int8_t>(Quadrant::NEL)] =
        new OctTree<POINT, MAX_POINTS>(p0, p1, _eps);

    // create child NWL
    p0[0] = _ll[0];
    p1[0] = x_mid;
    _children[static_cast<std::int8_t>(Quadrant::NWL)] =
        new OctTree<POINT, MAX_POINTS>(p0, p1, _eps);

    // create child SWL
    p0[1] = _ll[1];
    p1[1] = y_mid;
    _children[static_cast<std::int8_t>(Quadrant::SWL)] =
        new OctTree<POINT, MAX_POINTS>(_ll, p1, _eps);

    // create child NEU
    _children[static_cast<std::int8_t>(Quadrant::NEU)] =
        new OctTree<POINT, MAX_POINTS>(p1, _ur, _eps);

    // create child SEL
    p0[0] = x_mid;
    p1[0] = _ur[0];
    _children[static_cast<std::int8_t>(Quadrant::SEL)] =
        new OctTree<POINT, MAX_POINTS>(p0, p1, _eps);

    // create child NWU
    p0[0] = _ll[0];
    p0[1] = y_mid;
    p0[2] = z_mid;
    p1[0] = x_mid;
    p1[1] = _ur[1];
    p1[2] = _ur[2];
    _children[static_cast<std::int8_t>(Quadrant::NWU)] =
        new OctTree<POINT, MAX_POINTS>(p0, p1, _eps);

    // create child SWU
    p0[1] = _ll[1];
    p1[1] = y_mid;
    _children[static_cast<std::int8_t>(Quadrant::SWU)] =
        new OctTree<POINT, MAX_POINTS>(p0, p1, _eps);

    // create child SEU
    p0[0] = x_mid;
    p1[0] = _ur[0];
    p1[1] = y_mid;
    p1[2] = _ur[2];
    _children[static_cast<std::int8_t>(Quadrant::SEU)] =
        new OctTree<POINT, MAX_POINTS>(p0, p1, _eps);

    // add the passed point pnt to the children at first
    for (std::size_t k(0); k < 8; k++)
    {
        if (_children[k]->addPointToChild(pnt))
        {
            break;
        }
    }

    // distribute points to sub quadtrees
    const std::size_t n_pnts(_pnts.size());
    for (std::size_t j(0); j < n_pnts; j++)
    {
        for (auto c : _children)
        {
            if (c->addPointToChild(_pnts[j]))
            {
                break;
            }
        }
    }
    _is_leaf = false;
}

template <typename POINT, std::size_t MAX_POINTS>
bool OctTree<POINT, MAX_POINTS>::isOutside(POINT* pnt) const
{
    if ((*pnt)[0] < _ll[0] || (*pnt)[1] < _ll[1] || (*pnt)[2] < _ll[2])
    {
        return true;
    }
    if ((*pnt)[0] >= _ur[0] || (*pnt)[1] >= _ur[1] || (*pnt)[2] >= _ur[2])
    {
        return true;
    }
    return false;
}
}  // end namespace GeoLib
