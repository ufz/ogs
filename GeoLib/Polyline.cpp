/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-21
 * \brief  Implementation of the Polyline class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Polyline.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include "BaseLib/Error.h"
#include "AnalyticalGeometry.h"
#include "MathLib/GeometricBasics.h"

namespace GeoLib
{
Polyline::Polyline(const std::vector<Point*>& pnt_vec) :
    GeoObject(), _ply_pnts(pnt_vec)
{
    _length.push_back (0.0);
}

Polyline::Polyline(const Polyline& ply)
    : GeoObject(),
      _ply_pnts(ply._ply_pnts),
      _ply_pnt_ids(ply._ply_pnt_ids),
      _length(ply._length)
{}

void Polyline::write(std::ostream &os) const
{
    std::size_t size(_ply_pnt_ids.size());
    for (std::size_t k(0); k < size; k++)
        os << *(_ply_pnts[_ply_pnt_ids[k]]) << "\n";
}

bool Polyline::addPoint(std::size_t pnt_id)
{
    assert(pnt_id < _ply_pnts.size());
    std::size_t const n_pnts(_ply_pnt_ids.size());

    // don't insert point if this would result in identical IDs for two adjacent points
    if (n_pnts > 0 && _ply_pnt_ids.back() == pnt_id)
        return false;

    _ply_pnt_ids.push_back(pnt_id);

    if (n_pnts > 0)
    {
        double const act_dist(std::sqrt(MathLib::sqrDist(
            *_ply_pnts[_ply_pnt_ids[n_pnts-1]], *_ply_pnts[pnt_id])));
        double dist_until_now(0.0);
        if (n_pnts > 1)
            dist_until_now = _length[n_pnts - 1];

        _length.push_back(dist_until_now + act_dist);
    }
    return true;
}

bool Polyline::insertPoint(std::size_t pos, std::size_t pnt_id)
{
    assert(pnt_id < _ply_pnts.size());
    assert(pos <= _ply_pnt_ids.size());

    if (pos == _ply_pnt_ids.size()) {
        return addPoint(pnt_id);
    }

    // check if inserting pnt_id would result in two identical IDs for adjacent points
    if (pos == 0 && pnt_id == _ply_pnt_ids[0]) {
        return false;
    }
    if (pos != 0)
    {
        if (pos == (_ply_pnt_ids.size() - 1) && pnt_id == _ply_pnt_ids[pos])
        {
            return false;
        }
        if (pnt_id == _ply_pnt_ids[pos - 1] || pnt_id == _ply_pnt_ids[pos])
        {
            return false;
        }
        }

    auto const pos_dt(
        static_cast<std::vector<std::size_t>::difference_type>(pos));
    auto it(_ply_pnt_ids.begin() + pos_dt);
    _ply_pnt_ids.insert(it, pnt_id);

    if (_ply_pnt_ids.size() > 1) {
        // update the _length vector
        if (pos == 0) {
            // insert at first position
            double const act_dist(std::sqrt(MathLib::sqrDist(
                *_ply_pnts[_ply_pnt_ids[1]], *_ply_pnts[pnt_id])));
            _length.insert(_length.begin() + 1, act_dist);
            const std::size_t s(_length.size());
            for (std::size_t k(2); k < s; k++)
                _length[k] += _length[1];
        } else {
            if (pos == _ply_pnt_ids.size() - 1) {
                // insert at last position
                double const act_dist(std::sqrt(MathLib::sqrDist(
                    *_ply_pnts[_ply_pnt_ids[_ply_pnt_ids.size() - 2]],
                    *_ply_pnts[pnt_id])));
                double dist_until_now (0.0);
                if (_ply_pnt_ids.size() > 2)
                    dist_until_now = _length[_ply_pnt_ids.size() - 2];

                _length.insert(_length.begin() + pos_dt,
                               dist_until_now + act_dist);
            } else {
                // insert at arbitrary position within the vector
                double dist_until_now (0.0);
                if (pos > 1)
                    dist_until_now = _length[pos - 1];
                double len_seg0(std::sqrt(MathLib::sqrDist(
                                             *_ply_pnts[_ply_pnt_ids[pos - 1]],
                                             *_ply_pnts[pnt_id])));
                double len_seg1(std::sqrt(MathLib::sqrDist(
                                             *_ply_pnts[_ply_pnt_ids[pos + 1]],
                                             *_ply_pnts[pnt_id])));
                double update_dist(
                        len_seg0 + len_seg1 - (_length[pos] - dist_until_now));
                _length[pos] = dist_until_now + len_seg0;
                auto it1(_length.begin() + pos_dt + 1);
                _length.insert(it1, _length[pos] + len_seg1);
                for (it1 = _length.begin() + pos_dt + 2; it1 != _length.end();
                     ++it1)
                    *it1 += update_dist;
            }
        }
    }
    return true;
}

void Polyline::removePoint(std::size_t pos)
{
    if (pos >= _ply_pnt_ids.size())
        return;

    auto const pos_dt(
        static_cast<std::vector<std::size_t>::difference_type>(pos));
    _ply_pnt_ids.erase(_ply_pnt_ids.begin() + pos_dt);

    if (pos == _ply_pnt_ids.size())
    {
        _length.erase(_length.begin() + pos_dt);
        return;
    }

    const std::size_t n_ply_pnt_ids(_ply_pnt_ids.size());
    if (pos == 0) {
        double seg_length(_length[0]);
        for (std::size_t k(0); k < n_ply_pnt_ids; k++)
            _length[k] = _length[k + 1] - seg_length;
        _length.pop_back();
    } else {
        const double len_seg0(_length[pos] - _length[pos - 1]);
        const double len_seg1(_length[pos + 1] - _length[pos]);
        _length.erase(_length.begin() + pos_dt);
        const double len_new_seg(std::sqrt(MathLib::sqrDist(*_ply_pnts[_ply_pnt_ids[pos - 1]],
                                                       *_ply_pnts[_ply_pnt_ids[pos]])));
        double seg_length_diff(len_new_seg - len_seg0 - len_seg1);

        for (std::size_t k(pos); k < n_ply_pnt_ids; k++)
            _length[k] += seg_length_diff;
    }
}

std::size_t Polyline::getNumberOfPoints() const
{
    return _ply_pnt_ids.size();
}

std::size_t Polyline::getNumberOfSegments() const
{
    return _ply_pnt_ids.empty() ? 0 : _ply_pnt_ids.size()-1;
}

bool Polyline::isClosed() const
{
    if (_ply_pnt_ids.size() < 3)
        return false;

    if (_ply_pnt_ids.front() == _ply_pnt_ids.back())
        return true;

    return false;
}

bool Polyline::isCoplanar() const
{
    std::size_t const n_points (_ply_pnt_ids.size());
    if (n_points < 4)
        return true;

    GeoLib::Point const& p0 (*this->getPoint(0));
    GeoLib::Point const& p1 (*this->getPoint(1));
    GeoLib::Point const& p2 (*this->getPoint(2));
    for (std::size_t i=3; i<n_points; ++i)
    {
        if (!MathLib::isCoplanar(p0, p1, p2, *this->getPoint(i)))
        {
            DBUG ("Point %d is not coplanar to the first three points of the line.", i);
            return false;
        }
    }
    return true;
}

bool Polyline::isPointIDInPolyline(std::size_t pnt_id) const
{
    return std::find(_ply_pnt_ids.begin(), _ply_pnt_ids.end(), pnt_id) != _ply_pnt_ids.end();
}

std::size_t Polyline::getPointID(std::size_t i) const
{
    assert(i < _ply_pnt_ids.size());
    return _ply_pnt_ids[i];
}

LineSegment const Polyline::getSegment(std::size_t i) const
{
    assert(i < getNumberOfSegments());
    return LineSegment(_ply_pnts[_ply_pnt_ids[i]],
                       _ply_pnts[_ply_pnt_ids[i + 1]], false);
}

LineSegment Polyline::getSegment(std::size_t i)
{
    assert(i < getNumberOfSegments());
    return LineSegment(_ply_pnts[_ply_pnt_ids[i]],
                       _ply_pnts[_ply_pnt_ids[i + 1]], false);
}

void Polyline::setPointID(std::size_t idx, std::size_t id)
{
    assert(idx < _ply_pnt_ids.size());
    _ply_pnt_ids[idx] = id;
}

const Point* Polyline::getPoint(std::size_t i) const
{
    assert(i < _ply_pnt_ids.size());
    return _ply_pnts[_ply_pnt_ids[i]];
}

std::vector<Point*> const& Polyline::getPointsVec () const
{
    return _ply_pnts;
}

double Polyline::getLength (std::size_t k) const
{
    assert(k < _length.size());
    return _length[k];
}

Polyline* Polyline::constructPolylineFromSegments(const std::vector<Polyline*> &ply_vec,
                                                  double prox)
{
    std::size_t nLines = ply_vec.size();

    auto* new_ply = new Polyline(*ply_vec[0]);
    std::vector<GeoLib::Point*> pnt_vec(new_ply->getPointsVec());

    std::vector<Polyline*> local_ply_vec;
    for (std::size_t i = 1; i < nLines; i++)
        local_ply_vec.push_back(ply_vec[i]);

    while (!local_ply_vec.empty())
    {
        bool ply_found(false);
        prox *= prox; // square distance once to save time later
        for (auto it = local_ply_vec.begin(); it != local_ply_vec.end(); ++it)
        {
            if (pnt_vec == (*it)->getPointsVec())
            {
                std::size_t nPoints((*it)->getNumberOfPoints());

                //if (new_ply->getPointID(0) == (*it)->getPointID(0))
                if (pointsAreIdentical(pnt_vec, new_ply->getPointID(0),
                                       (*it)->getPointID(0), prox))
                {
                    auto* tmp = new Polyline((*it)->getPointsVec());
                    for (std::size_t k = 0; k < nPoints; k++)
                        tmp->addPoint((*it)->getPointID(nPoints - k - 1));

                    std::size_t new_ply_size(new_ply->getNumberOfPoints());
                    for (std::size_t k = 1; k < new_ply_size; k++)
                        tmp->addPoint(new_ply->getPointID(k));
                    delete new_ply;
                    new_ply = tmp;
                    ply_found = true;
                }
                //else if (new_ply->getPointID(0) == (*it)->getPointID(nPoints-1))
                else if (pointsAreIdentical(pnt_vec, new_ply->getPointID(0),
                                            (*it)->getPointID(nPoints - 1), prox))
                {
                    auto* tmp = new Polyline(**it);
                    std::size_t new_ply_size(new_ply->getNumberOfPoints());
                    for (std::size_t k = 1; k < new_ply_size; k++)
                        tmp->addPoint(new_ply->getPointID(k));
                    delete new_ply;
                    new_ply = tmp;
                    ply_found = true;
                }
                //else if (new_ply->getPointID(new_ply->getNumberOfPoints()-1) == (*it)->getPointID(0))
                else if (pointsAreIdentical(pnt_vec,
                                            new_ply->getPointID(new_ply->
                                                                getNumberOfPoints()
                                                                - 1),
                                            (*it)->getPointID(0), prox))
                {
                    for (std::size_t k = 1; k < nPoints; k++)
                        new_ply->addPoint((*it)->getPointID(k));
                    ply_found = true;
                }
                //else if (new_ply->getPointID(new_ply->getNumberOfPoints()-1) == (*it)->getPointID(nPoints-1))
                else if (pointsAreIdentical(pnt_vec,
                                            new_ply->getPointID(new_ply->
                                                                getNumberOfPoints()
                                                                - 1),
                                            (*it)->getPointID(nPoints - 1), prox))
                {
                    for (std::size_t k = 1; k < nPoints; k++)
                        new_ply->addPoint((*it)->getPointID(nPoints - k - 1));
                    ply_found = true;
                }
                if (ply_found)
                {
                    local_ply_vec.erase(it);
                    break;
                }
            }
            else
                ERR("Error in Polyline::contructPolylineFromSegments() - Line segments use different point vectors.");
        }

        if (!ply_found)
        {
            ERR("Error in Polyline::contructPolylineFromSegments() - Not all segments are connected.");
            new_ply = nullptr;
            break;
        }
    }
    return new_ply;
}

void Polyline::closePolyline()
{
    if (getNumberOfPoints() < 2) {
        ERR("Polyline::closePolyline(): Input polyline needs to be composed of at least three points.");
    }
    if (! isClosed()) {
        addPoint(getPointID(0));
    }
}

Location Polyline::getLocationOfPoint (std::size_t k, GeoLib::Point const & pnt) const
{
    assert (k < _ply_pnt_ids.size() - 1);

    GeoLib::Point const& source (*(_ply_pnts[_ply_pnt_ids[k]]));
    GeoLib::Point const& dest (*(_ply_pnts[_ply_pnt_ids[k + 1]]));
    long double a[2] = {dest[0] - source[0], dest[1] - source[1]}; // vector
    long double b[2] = {pnt[0] - source[0], pnt[1] - source[1]}; // vector

    long double det_2x2 (a[0] * b[1] - a[1] * b[0]);

    if (det_2x2 > std::numeric_limits<double>::epsilon())
        return Location::LEFT;
    if (std::numeric_limits<double>::epsilon() < std::abs(det_2x2))
        return Location::RIGHT;
    if (a[0] * b[0] < 0.0 || a[1] * b[1] < 0.0)
        return Location::BEHIND;
    if (a[0] * a[0] + a[1] * a[1] < b[0] * b[0] + b[1] * b[1])
        return Location::BEYOND;
    if (MathLib::sqrDist (pnt,
                          *_ply_pnts[_ply_pnt_ids[k]]) < pow(std::numeric_limits<double>::epsilon(),2))
        return Location::SOURCE;
    if (MathLib::sqrDist (pnt,
                          *_ply_pnts[_ply_pnt_ids[k + 1]]) <
        std::sqrt(std::numeric_limits<double>::epsilon()))
        return Location::DESTINATION;
    return Location::BETWEEN;
}

void Polyline::updatePointIDs(const std::vector<std::size_t> &pnt_ids)
{
    for (auto it = this->_ply_pnt_ids.begin(); it!=this->_ply_pnt_ids.end();)
    {
        if (pnt_ids[*it] != *it)
        {
            if (it!=this->_ply_pnt_ids.begin() && (pnt_ids[*it] == pnt_ids[*(it-1)]))
                it = this->_ply_pnt_ids.erase(it);
            else
            {
                *it = pnt_ids[*it];
                ++it;
            }
        }
        else
            ++it;
    }
}

double Polyline::getDistanceAlongPolyline(const MathLib::Point3d& pnt,
    const double epsilon_radius) const
{
    double dist(-1.0), lambda;
    bool found = false;
    // loop over all line segments of the polyline
    for (std::size_t k = 0; k < getNumberOfSegments(); k++) {
        // is the orthogonal projection of the j-th node to the
        // line g(lambda) = _ply->getPoint(k) + lambda * (_ply->getPoint(k+1) - _ply->getPoint(k))
        // at the k-th line segment of the polyline, i.e. 0 <= lambda <= 1?
        if (MathLib::calcProjPntToLineAndDists(pnt.getCoords(),
                        (getPoint(k))->getCoords(), (getPoint(k + 1))->getCoords(),
                        lambda, dist) <= epsilon_radius) {

            double act_length_of_ply(getLength(k));
            double seg_length (getLength(k+1)-getLength(k));
            double lower_lambda (- epsilon_radius / seg_length);
            double upper_lambda (1 + epsilon_radius / seg_length);

            if (lower_lambda <= lambda && lambda <= upper_lambda) {
                found = true;
                dist = act_length_of_ply + dist;
                break;
            } // end if lambda
        }
    } // end line segment loop

    if (! found)
        dist = -1.0;
    return dist;
}

Polyline::SegmentIterator::SegmentIterator(Polyline const& polyline,
                                           std::size_t segment_number)
    : _polyline(&polyline),
      _segment_number(
          static_cast<std::vector<GeoLib::Point*>::size_type>(segment_number))
{}

Polyline::SegmentIterator::SegmentIterator(SegmentIterator const& src)
    : _polyline(src._polyline), _segment_number(src._segment_number)
{}

Polyline::SegmentIterator& Polyline::SegmentIterator::operator=(
    SegmentIterator const& rhs)
{
    if (&rhs == this)
        return *this;

    _polyline = rhs._polyline;
    _segment_number = rhs._segment_number;
    return *this;
}

std::size_t Polyline::SegmentIterator::getSegmentNumber() const
{
    return static_cast<std::size_t>(_segment_number);
}

Polyline::SegmentIterator& Polyline::SegmentIterator::operator++()
{
    ++_segment_number;
    return *this;
}

LineSegment const Polyline::SegmentIterator::operator*() const
{
    return _polyline->getSegment(_segment_number);
}

LineSegment Polyline::SegmentIterator::operator*()
{
    return _polyline->getSegment(_segment_number);
}

bool Polyline::SegmentIterator::operator==(SegmentIterator const& other)
{
    return !(*this != other);
}

bool Polyline::SegmentIterator::operator!=(SegmentIterator const& other)
{
    return other._segment_number != _segment_number ||
           other._polyline != _polyline;
}

Polyline::SegmentIterator& Polyline::SegmentIterator::operator+=(
    std::vector<GeoLib::Point>::difference_type n)
{
    if (n < 0) {
        _segment_number -=
            static_cast<std::vector<GeoLib::Point>::size_type>(-n);
    } else {
        _segment_number +=
            static_cast<std::vector<GeoLib::Point>::size_type>(n);
    }
    if (_segment_number > _polyline->getNumberOfSegments())
        OGS_FATAL("");
    return *this;
}

Polyline::SegmentIterator Polyline::SegmentIterator::operator+(
    std::vector<GeoLib::Point>::difference_type n)
{
    SegmentIterator t(*this);
    t += n;
    return t;
}

Polyline::SegmentIterator& Polyline::SegmentIterator::operator-=(
    std::vector<GeoLib::Point>::difference_type n)
{
    if (n >= 0) {
        _segment_number -=
            static_cast<std::vector<GeoLib::Point>::size_type>(n);
    } else {
        _segment_number +=
            static_cast<std::vector<GeoLib::Point>::size_type>(-n);
    }
    if (_segment_number > _polyline->getNumberOfSegments())
        OGS_FATAL("");
    return *this;
}

Polyline::SegmentIterator Polyline::SegmentIterator::operator-(
    std::vector<GeoLib::Point>::difference_type n)
{
    Polyline::SegmentIterator t(*this);
    t -= n;
    return t;
}

std::ostream& operator<< (std::ostream &os, const Polyline &pl)
{
    pl.write (os);
    return os;
}

bool containsEdge (const Polyline& ply, std::size_t id0, std::size_t id1)
{
    if (id0 == id1)
    {
        ERR("no valid edge id0 == id1 == %d.", id0);
        return false;
    }
    if (id0 > id1)
        std::swap (id0,id1);
    const std::size_t n (ply.getNumberOfPoints() - 1);
    for (std::size_t k(0); k < n; k++)
    {
        std::size_t ply_pnt0 (ply.getPointID (k));
        std::size_t ply_pnt1 (ply.getPointID (k + 1));
        if (ply_pnt0 > ply_pnt1)
            std::swap (ply_pnt0, ply_pnt1);
        if (ply_pnt0 == id0 && ply_pnt1 == id1)
            return true;
    }
    return false;
}

bool operator==(Polyline const& lhs, Polyline const& rhs)
{
    if (lhs.getNumberOfPoints() != rhs.getNumberOfPoints())
        return false;

    const std::size_t n(lhs.getNumberOfPoints());
    for (std::size_t k(0); k < n; k++)
        if (lhs.getPointID(k) != rhs.getPointID(k))
            return false;

    return true;
}

bool pointsAreIdentical(const std::vector<Point*> &pnt_vec,
                        std::size_t i,
                        std::size_t j,
                        double prox)
{
    if (i == j)
        return true;
    return MathLib::sqrDist(*pnt_vec[i], *pnt_vec[j]) < prox;
}
} // end namespace GeoLib
