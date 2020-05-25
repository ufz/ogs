/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-21
 * \brief  Implementation of the Polyline class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Polyline.h"

#include <algorithm>
#include "AnalyticalGeometry.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "MathLib/GeometricBasics.h"

namespace GeoLib
{
Polyline::Polyline(const std::vector<Point*>& pnt_vec) : ply_pnts_(pnt_vec)
{
    length_.push_back (0.0);
}

Polyline::Polyline(const Polyline& ply)
    : ply_pnts_(ply.ply_pnts_),
      ply_pnt_ids_(ply.ply_pnt_ids_),
      length_(ply.length_)
{}

void Polyline::write(std::ostream &os) const
{
    std::size_t size(ply_pnt_ids_.size());
    for (std::size_t k(0); k < size; k++)
    {
        os << *(ply_pnts_[ply_pnt_ids_[k]]) << "\n";
    }
}

bool Polyline::addPoint(std::size_t pnt_id)
{
    assert(pnt_id < ply_pnts_.size());
    std::size_t const n_pnts(ply_pnt_ids_.size());

    // don't insert point if this would result in identical IDs for two adjacent points
    if (n_pnts > 0 && ply_pnt_ids_.back() == pnt_id)
    {
        return false;
    }

    ply_pnt_ids_.push_back(pnt_id);

    if (n_pnts > 0)
    {
        double const act_dist(std::sqrt(MathLib::sqrDist(
            *ply_pnts_[ply_pnt_ids_[n_pnts-1]], *ply_pnts_[pnt_id])));
        double dist_until_now(0.0);
        if (n_pnts > 1)
        {
            dist_until_now = length_[n_pnts - 1];
        }

        length_.push_back(dist_until_now + act_dist);
    }
    return true;
}

bool Polyline::insertPoint(std::size_t pos, std::size_t pnt_id)
{
    assert(pnt_id < ply_pnts_.size());
    assert(pos <= ply_pnt_ids_.size());

    if (pos == ply_pnt_ids_.size()) {
        return addPoint(pnt_id);
    }

    // check if inserting pnt_id would result in two identical IDs for adjacent points
    if (pos == 0 && pnt_id == ply_pnt_ids_[0]) {
        return false;
    }
    if (pos != 0)
    {
        if (pos == (ply_pnt_ids_.size() - 1) && pnt_id == ply_pnt_ids_[pos])
        {
            return false;
        }
        if (pnt_id == ply_pnt_ids_[pos - 1] || pnt_id == ply_pnt_ids_[pos])
        {
            return false;
        }
        }

    auto const pos_dt(
        static_cast<std::vector<std::size_t>::difference_type>(pos));
    auto it(ply_pnt_ids_.begin() + pos_dt);
    ply_pnt_ids_.insert(it, pnt_id);

    if (ply_pnt_ids_.size() > 1) {
        // update the length_ vector
        if (pos == 0) {
            // insert at first position
            double const act_dist(std::sqrt(MathLib::sqrDist(
                *ply_pnts_[ply_pnt_ids_[1]], *ply_pnts_[pnt_id])));
            length_.insert(length_.begin() + 1, act_dist);
            const std::size_t s(length_.size());
            for (std::size_t k(2); k < s; k++)
            {
                length_[k] += length_[1];
            }
        } else {
            if (pos == ply_pnt_ids_.size() - 1) {
                // insert at last position
                double const act_dist(std::sqrt(MathLib::sqrDist(
                    *ply_pnts_[ply_pnt_ids_[ply_pnt_ids_.size() - 2]],
                    *ply_pnts_[pnt_id])));
                double dist_until_now (0.0);
                if (ply_pnt_ids_.size() > 2)
                {
                    dist_until_now = length_[ply_pnt_ids_.size() - 2];
                }

                length_.insert(length_.begin() + pos_dt,
                               dist_until_now + act_dist);
            } else {
                // insert at arbitrary position within the vector
                double dist_until_now (0.0);
                if (pos > 1)
                {
                    dist_until_now = length_[pos - 1];
                }
                double len_seg0(std::sqrt(MathLib::sqrDist(
                                             *ply_pnts_[ply_pnt_ids_[pos - 1]],
                                             *ply_pnts_[pnt_id])));
                double len_seg1(std::sqrt(MathLib::sqrDist(
                                             *ply_pnts_[ply_pnt_ids_[pos + 1]],
                                             *ply_pnts_[pnt_id])));
                double update_dist(
                        len_seg0 + len_seg1 - (length_[pos] - dist_until_now));
                length_[pos] = dist_until_now + len_seg0;
                auto it1(length_.begin() + pos_dt + 1);
                length_.insert(it1, length_[pos] + len_seg1);
                for (it1 = length_.begin() + pos_dt + 2; it1 != length_.end();
                     ++it1)
                {
                    *it1 += update_dist;
                }
            }
        }
    }
    return true;
}

void Polyline::removePoint(std::size_t pos)
{
    if (pos >= ply_pnt_ids_.size())
    {
        return;
    }

    auto const pos_dt(
        static_cast<std::vector<std::size_t>::difference_type>(pos));
    ply_pnt_ids_.erase(ply_pnt_ids_.begin() + pos_dt);

    if (pos == ply_pnt_ids_.size())
    {
        length_.erase(length_.begin() + pos_dt);
        return;
    }

    const std::size_t n_ply_pnt_ids(ply_pnt_ids_.size());
    if (pos == 0) {
        double seg_length(length_[0]);
        for (std::size_t k(0); k < n_ply_pnt_ids; k++)
        {
            length_[k] = length_[k + 1] - seg_length;
        }
        length_.pop_back();
    } else {
        const double len_seg0(length_[pos] - length_[pos - 1]);
        const double len_seg1(length_[pos + 1] - length_[pos]);
        length_.erase(length_.begin() + pos_dt);
        const double len_new_seg(std::sqrt(MathLib::sqrDist(*ply_pnts_[ply_pnt_ids_[pos - 1]],
                                                       *ply_pnts_[ply_pnt_ids_[pos]])));
        double seg_length_diff(len_new_seg - len_seg0 - len_seg1);

        for (std::size_t k(pos); k < n_ply_pnt_ids; k++)
        {
            length_[k] += seg_length_diff;
        }
    }
}

std::size_t Polyline::getNumberOfPoints() const
{
    return ply_pnt_ids_.size();
}

std::size_t Polyline::getNumberOfSegments() const
{
    return ply_pnt_ids_.empty() ? 0 : ply_pnt_ids_.size()-1;
}

bool Polyline::isClosed() const
{
    if (ply_pnt_ids_.size() < 3)
    {
        return false;
    }

    return ply_pnt_ids_.front() == ply_pnt_ids_.back();
}

bool Polyline::isCoplanar() const
{
    std::size_t const n_points (ply_pnt_ids_.size());
    if (n_points < 4)
    {
        return true;
    }

    GeoLib::Point const& p0 (*this->getPoint(0));
    GeoLib::Point const& p1 (*this->getPoint(1));
    GeoLib::Point const& p2 (*this->getPoint(2));
    for (std::size_t i=3; i<n_points; ++i)
    {
        if (!MathLib::isCoplanar(p0, p1, p2, *this->getPoint(i)))
        {
            DBUG(
                "Point {:d} is not coplanar to the first three points of the "
                "line.",
                i);
            return false;
        }
    }
    return true;
}

bool Polyline::isPointIDInPolyline(std::size_t pnt_id) const
{
    return std::find(ply_pnt_ids_.begin(), ply_pnt_ids_.end(), pnt_id) != ply_pnt_ids_.end();
}

std::size_t Polyline::getPointID(std::size_t i) const
{
    assert(i < ply_pnt_ids_.size());
    return ply_pnt_ids_[i];
}

LineSegment Polyline::getSegment(std::size_t i) const
{
    assert(i < getNumberOfSegments());
    return LineSegment(ply_pnts_[ply_pnt_ids_[i]],
                       ply_pnts_[ply_pnt_ids_[i + 1]], false);
}

LineSegment Polyline::getSegment(std::size_t i)
{
    assert(i < getNumberOfSegments());
    return LineSegment(ply_pnts_[ply_pnt_ids_[i]],
                       ply_pnts_[ply_pnt_ids_[i + 1]], false);
}

void Polyline::setPointID(std::size_t idx, std::size_t id)
{
    assert(idx < ply_pnt_ids_.size());
    ply_pnt_ids_[idx] = id;
}

const Point* Polyline::getPoint(std::size_t i) const
{
    assert(i < ply_pnt_ids_.size());
    return ply_pnts_[ply_pnt_ids_[i]];
}

std::vector<Point*> const& Polyline::getPointsVec () const
{
    return ply_pnts_;
}

double Polyline::getLength (std::size_t k) const
{
    assert(k < length_.size());
    return length_[k];
}

Polyline* Polyline::constructPolylineFromSegments(const std::vector<Polyline*> &ply_vec,
                                                  double prox)
{
    std::size_t nLines = ply_vec.size();

    auto* new_ply = new Polyline(*ply_vec[0]);
    std::vector<GeoLib::Point*> pnt_vec(new_ply->getPointsVec());

    std::vector<Polyline*> local_ply_vec;
    for (std::size_t i = 1; i < nLines; i++)
    {
        local_ply_vec.push_back(ply_vec[i]);
    }

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
                    {
                        tmp->addPoint((*it)->getPointID(nPoints - k - 1));
                    }

                    std::size_t new_ply_size(new_ply->getNumberOfPoints());
                    for (std::size_t k = 1; k < new_ply_size; k++)
                    {
                        tmp->addPoint(new_ply->getPointID(k));
                    }
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
                    {
                        tmp->addPoint(new_ply->getPointID(k));
                    }
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
                    {
                        new_ply->addPoint((*it)->getPointID(k));
                    }
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
                    {
                        new_ply->addPoint((*it)->getPointID(nPoints - k - 1));
                    }
                    ply_found = true;
                }
                if (ply_found)
                {
                    local_ply_vec.erase(it);
                    break;
                }
            }
            else
            {
                ERR("Error in Polyline::contructPolylineFromSegments() - Line segments use different point vectors.");
            }
        }

        if (!ply_found)
        {
            ERR("Error in Polyline::contructPolylineFromSegments() - Not all segments are connected.");
            delete new_ply;
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
    assert (k < ply_pnt_ids_.size() - 1);

    GeoLib::Point const& source (*(ply_pnts_[ply_pnt_ids_[k]]));
    GeoLib::Point const& dest (*(ply_pnts_[ply_pnt_ids_[k + 1]]));
    long double a[2] = {dest[0] - source[0], dest[1] - source[1]}; // vector
    long double b[2] = {pnt[0] - source[0], pnt[1] - source[1]}; // vector

    long double det_2x2 (a[0] * b[1] - a[1] * b[0]);

    if (det_2x2 > std::numeric_limits<double>::epsilon())
    {
        return Location::LEFT;
    }
    if (std::numeric_limits<double>::epsilon() < std::abs(det_2x2))
    {
        return Location::RIGHT;
    }
    if (a[0] * b[0] < 0.0 || a[1] * b[1] < 0.0)
    {
        return Location::BEHIND;
    }
    if (a[0] * a[0] + a[1] * a[1] < b[0] * b[0] + b[1] * b[1])
    {
        return Location::BEYOND;
    }
    if (MathLib::sqrDist(pnt, *ply_pnts_[ply_pnt_ids_[k]]) <
        pow(std::numeric_limits<double>::epsilon(), 2))
    {
        return Location::SOURCE;
    }
    if (MathLib::sqrDist(pnt, *ply_pnts_[ply_pnt_ids_[k + 1]]) <
        std::sqrt(std::numeric_limits<double>::epsilon()))
    {
        return Location::DESTINATION;
    }
    return Location::BETWEEN;
}

double Polyline::getDistanceAlongPolyline(const MathLib::Point3d& pnt,
    const double epsilon_radius) const
{
    double dist(-1.0);
    double lambda;
    bool found = false;
    // loop over all line segments of the polyline
    for (std::size_t k = 0; k < getNumberOfSegments(); k++) {
        // is the orthogonal projection of the j-th node to the
        // line g(lambda) = ply_->getPoint(k) + lambda * (ply_->getPoint(k+1) - ply_->getPoint(k))
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

    if (!found)
    {
        dist = -1.0;
    }
    return dist;
}

Polyline::SegmentIterator::SegmentIterator(Polyline const& polyline,
                                           std::size_t segment_number)
    : polyline_(&polyline),
      segment_number_(
          static_cast<std::vector<GeoLib::Point*>::size_type>(segment_number))
{}

Polyline::SegmentIterator::SegmentIterator(SegmentIterator const& src)
    : polyline_(src.polyline_), segment_number_(src.segment_number_)
{}

Polyline::SegmentIterator& Polyline::SegmentIterator::operator=(
    SegmentIterator const& rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    polyline_ = rhs.polyline_;
    segment_number_ = rhs.segment_number_;
    return *this;
}

std::size_t Polyline::SegmentIterator::getSegmentNumber() const
{
    return static_cast<std::size_t>(segment_number_);
}

Polyline::SegmentIterator& Polyline::SegmentIterator::operator++()
{
    ++segment_number_;
    return *this;
}

LineSegment Polyline::SegmentIterator::operator*() const
{
    return polyline_->getSegment(segment_number_);
}

LineSegment Polyline::SegmentIterator::operator*()
{
    return polyline_->getSegment(segment_number_);
}

bool Polyline::SegmentIterator::operator==(SegmentIterator const& other)
{
    return !(*this != other);
}

bool Polyline::SegmentIterator::operator!=(SegmentIterator const& other)
{
    return other.segment_number_ != segment_number_ ||
           other.polyline_ != polyline_;
}

Polyline::SegmentIterator& Polyline::SegmentIterator::operator+=(
    std::vector<GeoLib::Point>::difference_type n)
{
    if (n < 0) {
        segment_number_ -=
            static_cast<std::vector<GeoLib::Point>::size_type>(-n);
    } else {
        segment_number_ +=
            static_cast<std::vector<GeoLib::Point>::size_type>(n);
    }
    if (segment_number_ > polyline_->getNumberOfSegments())
    {
        OGS_FATAL("");
    }
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
        segment_number_ -=
            static_cast<std::vector<GeoLib::Point>::size_type>(n);
    } else {
        segment_number_ +=
            static_cast<std::vector<GeoLib::Point>::size_type>(-n);
    }
    if (segment_number_ > polyline_->getNumberOfSegments())
    {
        OGS_FATAL("");
    }
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
        ERR("no valid edge id0 == id1 == {:d}.", id0);
        return false;
    }
    if (id0 > id1)
    {
        std::swap(id0, id1);
    }
    const std::size_t n (ply.getNumberOfPoints() - 1);
    for (std::size_t k(0); k < n; k++)
    {
        std::size_t ply_pnt0 (ply.getPointID (k));
        std::size_t ply_pnt1 (ply.getPointID (k + 1));
        if (ply_pnt0 > ply_pnt1)
        {
            std::swap(ply_pnt0, ply_pnt1);
        }
        if (ply_pnt0 == id0 && ply_pnt1 == id1)
        {
            return true;
        }
    }
    return false;
}

bool operator==(Polyline const& lhs, Polyline const& rhs)
{
    if (lhs.getNumberOfPoints() != rhs.getNumberOfPoints())
    {
        return false;
    }

    const std::size_t n(lhs.getNumberOfPoints());
    for (std::size_t k(0); k < n; k++)
    {
        if (lhs.getPointID(k) != rhs.getPointID(k))
        {
            return false;
        }
    }

    return true;
}

bool pointsAreIdentical(const std::vector<Point*> &pnt_vec,
                        std::size_t i,
                        std::size_t j,
                        double prox)
{
    if (i == j)
    {
        return true;
    }
    return MathLib::sqrDist(*pnt_vec[i], *pnt_vec[j]) < prox;
}
} // end namespace GeoLib
