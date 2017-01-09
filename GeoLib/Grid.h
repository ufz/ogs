/**
 * \file
 * \author Thomas Fischer
 * \date   2012-02-02
 * \brief  Definition of the Grid class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GRID_H_
#define GRID_H_

#include <bitset>
#include <vector>

#include <logog/include/logog.hpp>

// GeoLib
#include "AABB.h"
#include "GEOObjects.h"

// MathLib
#include "MathLib/Point3d.h"
#include "MathLib/MathTools.h"

namespace GeoLib
{
template <typename POINT>
class Grid : public GeoLib::AABB
{
public:
    /**
     * @brief The constructor of the grid object takes a vector of points or nodes. Furthermore the
     * user can specify the *average* maximum number of points per grid cell.
     *
     * The number of grid cells are computed with the following formula
     * \f$\frac{n_{points}}{n_{cells}} \le n_{max\_per\_cell}\f$
     *
     * In order to limit the memory wasting the maximum number of points per grid cell
     * (in the average) should be a power of two (since std::vector objects resize itself
     * with this step size).
     *
     * @param first, last the range of elements to examine
     * @param items_per_cell (input) max number per grid cell in the average (default 512)
     *
     */
    template <typename InputIterator>
    Grid(InputIterator first, InputIterator last, std::size_t items_per_cell = 512);

    /**
     * This is the destructor of the class. It deletes the internal data structures *not*
     * including the pointers to the points.
     */
    virtual ~Grid()
    {
        delete [] _grid_cell_nodes_map;
    }

    /**
     * The method calculates the grid cell the given point is belonging to, i.e.,
     * the (internal) coordinates of the grid cell are computed. The method searches the actual
     * grid cell and all its neighbors for the POINT object which has the smallest
     * distance. A pointer to this object is returned.
     *
     * If there is not such a point, i.e., all the searched grid cells do not contain any
     * POINT object a nullptr is returned.
     *
     * @param pnt search point
     * @return a pointer to the point with the smallest distance within the grid cells that are
     * outlined above or nullptr
     */
    template <typename P> POINT* getNearestPoint(P const& pnt) const;

    template <typename P> std::vector<std::size_t> getPointsInEpsilonEnvironment(
        P const& pnt, double eps) const;

    /**
     * Method fetches the vectors of all grid cells intersecting the axis aligned cuboid
     * defined by two points. The first point with minimal coordinates in all directions.
     * The second point with maximal coordinates in all directions.
     *
     * @param center (input) the center point of the axis aligned cube
     * @param half_len (input) half of the edge length of the axis aligned cube
     * @return vector of vectors of points within grid cells that intersects
     * the axis aligned cube
     */
    template <typename P>
    std::vector<std::vector<POINT*> const*>
    getPntVecsOfGridCellsIntersectingCube(P const& center, double half_len) const;

    void getPntVecsOfGridCellsIntersectingCuboid(
        MathLib::Point3d const& min_pnt,
        MathLib::Point3d const& max_pnt,
        std::vector<std::vector<POINT*> const*>& pnts) const;

#ifndef NDEBUG
    /**
     * Method creates a geometry for every mesh grid box. Additionally it
     * creates one geometry containing all the box geometries.
     * @param geo_obj
     */
    void createGridGeometry(GeoLib::GEOObjects* geo_obj) const;
#endif

private:
    /// Computes the number of grid cells per spatial dimension the objects
    /// (points or mesh nodes) will be sorted in.
    /// On the one hand the number of grid cells should be small to reduce the
    /// management overhead. On the other hand the number should be large such
    /// that each grid cell contains only a small number of objects.
    /// Under the assumption that the points are distributed equidistant in
    /// space the grid cells should be as cubical as possible.
    /// At first it is necessary to determine the spatial dimension the grid
    /// should have. The dimensions are computed from the spatial extensions
    /// Let \f$\max\f$ be the largest spatial extension. The grid will have a
    /// spatial dimension if the ratio of the corresponding spatial extension
    /// and the maximal extension is \f$\ge 10^{-4}\f$.
    /// The second step consists of computing the number of cells per dimension.
    void initNumberOfSteps(std::size_t n_per_cell,
        std::size_t n_pnts, std::array<double,3> const& extensions);

    /**
     * Method calculates the grid cell coordinates for the given point pnt. If
     * the point is located outside of the bounding box the coordinates of the
     * grid cell on the border that is nearest to given point will be returned.
     * @param pnt (input) the coordinates of the point
     * @return the coordinates of the grid cell
     */
    template <typename T>
    std::array<std::size_t,3> getGridCoords(T const& pnt) const;

    /**
     *
     * point numbering of the grid cell is as follow
     * @code
     *         7 -------- 6
     *        /:         /|
     *       / :        / |
     *      /  :       /  |
     *     /   :      /   |
     *    4 -------- 5    |
     *    |    3 ....|... 2
     *    |   .      |   /
     *    |  .       |  /
     *    | .        | /
     *    |.         |/
     *    0 -------- 1
     *    @endcode
     *    the face numbering is as follow:
     *    face    nodes
     *    0        0,3,2,1 bottom
     *    1        0,1,5,4 front
     *    2        1,2,6,5 right
     *    3        2,3,7,6 back
     *    4        3,0,4,7 left
     *    5        4,5,6,7 top
     * @param pnt (input) coordinates of the point
     * @param coordinates of the grid cell
     * @return squared distances of the point to the faces of the grid cell
     * ordered in the same sequence as above described
     */
    template <typename P>
    std::array<double, 6> getPointCellBorderDistances(
        P const& pnt, std::array<std::size_t,3> const& coordinates) const;

    template <typename P>
    bool calcNearestPointInGridCell(P const& pnt,
        std::array<std::size_t,3> const& coords,
        double &sqr_min_dist,
        POINT* &nearest_pnt) const;

    static POINT* copyOrAddress(POINT& p) { return &p; }
    static POINT const* copyOrAddress(POINT const& p) { return &p; }
    static POINT* copyOrAddress(POINT* p) { return p; }

    std::array<std::size_t,3> _n_steps;
    std::array<double, 3> _step_sizes;
    std::array<double, 3> _inverse_step_sizes;
    /**
     * This is an array that stores pointers to POINT objects.
     */
    std::vector<POINT*>* _grid_cell_nodes_map;
};

template <typename POINT>
template <typename InputIterator>
Grid<POINT>::Grid(InputIterator first, InputIterator last,
    std::size_t max_num_per_grid_cell)
    : GeoLib::AABB(first, last), _n_steps({{1,1,1}}),
        _step_sizes({{0.0,0.0,0.0}}), _inverse_step_sizes({{0.0,0.0,0.0}}),
        _grid_cell_nodes_map(nullptr)
{
    auto const n_pnts(std::distance(first,last));

    std::array<double, 3> delta = {{_max_pnt[0] - _min_pnt[0],
                                    _max_pnt[1] - _min_pnt[1],
                                    _max_pnt[2] - _min_pnt[2]}};

    // enlarge delta
    for (auto & d : delta)
        d = std::nextafter(d, std::numeric_limits<double>::max());

    assert(n_pnts > 0);
    initNumberOfSteps(max_num_per_grid_cell, static_cast<std::size_t>(n_pnts), delta);

    const std::size_t n_plane(_n_steps[0] * _n_steps[1]);
    _grid_cell_nodes_map = new std::vector<POINT*>[n_plane * _n_steps[2]];

    // some frequently used expressions to fill the grid vectors
    for (std::size_t k(0); k < 3; k++) {
        if (std::abs(delta[k]) < std::numeric_limits<double>::epsilon()) {
            delta[k] = std::numeric_limits<double>::epsilon();
        }
        _step_sizes[k] = delta[k] / _n_steps[k];
        _inverse_step_sizes[k] = 1.0 / _step_sizes[k];
    }

    // fill the grid vectors
    InputIterator it(first);
    while (it != last) {
        std::array<std::size_t,3> coords(getGridCoords(*copyOrAddress(*it)));
        if (coords < _n_steps) {
            std::size_t const pos(coords[0]+coords[1]*_n_steps[0]+coords[2]*n_plane);
            _grid_cell_nodes_map[pos].push_back(const_cast<POINT*>(copyOrAddress(*it)));
        } else {
            ERR("Grid constructor: error computing indices [%d, %d, %d], "
                "max indices [%d, %d, %d].", coords[0], coords[1], coords[2],
                _n_steps[0], _n_steps[1], _n_steps[2]);
        }
        it++;
    }
}

template<typename POINT>
template <typename P>
std::vector<std::vector<POINT*> const*>
Grid<POINT>::getPntVecsOfGridCellsIntersectingCube(P const& center,
    double half_len) const
{
    std::vector<std::vector<POINT*> const*> pnts;
    MathLib::Point3d tmp_pnt{
        {{center[0]-half_len, center[1]-half_len, center[2]-half_len}}}; // min
    std::array<std::size_t,3> min_coords(getGridCoords(tmp_pnt));

    tmp_pnt[0] = center[0] + half_len;
    tmp_pnt[1] = center[1] + half_len;
    tmp_pnt[2] = center[2] + half_len;
    std::array<std::size_t,3> max_coords(getGridCoords(tmp_pnt));

    std::size_t coords[3], steps0_x_steps1(_n_steps[0] * _n_steps[1]);
    for (coords[0] = min_coords[0]; coords[0] < max_coords[0] + 1; coords[0]++) {
        for (coords[1] = min_coords[1]; coords[1] < max_coords[1] + 1; coords[1]++) {
            const std::size_t coords0_p_coords1_x_steps0(coords[0] + coords[1] * _n_steps[0]);
            for (coords[2] = min_coords[2]; coords[2] < max_coords[2] + 1; coords[2]++) {
                pnts.push_back(&(_grid_cell_nodes_map[coords0_p_coords1_x_steps0 + coords[2]
                                * steps0_x_steps1]));
            }
        }
    }
    return pnts;
}

template<typename POINT>
void Grid<POINT>::getPntVecsOfGridCellsIntersectingCuboid(
    MathLib::Point3d const& min_pnt,
    MathLib::Point3d const& max_pnt,
    std::vector<std::vector<POINT*> const*>& pnts) const
{
    std::array<std::size_t,3> min_coords(getGridCoords(min_pnt));
    std::array<std::size_t,3> max_coords(getGridCoords(max_pnt));

    std::size_t coords[3], steps0_x_steps1(_n_steps[0] * _n_steps[1]);
    for (coords[0] = min_coords[0]; coords[0] < max_coords[0] + 1; coords[0]++) {
        for (coords[1] = min_coords[1]; coords[1] < max_coords[1] + 1; coords[1]++) {
            const std::size_t coords0_p_coords1_x_steps0(coords[0] + coords[1] * _n_steps[0]);
            for (coords[2] = min_coords[2]; coords[2] < max_coords[2] + 1; coords[2]++) {
                pnts.push_back(&(_grid_cell_nodes_map[coords0_p_coords1_x_steps0 + coords[2]
                                * steps0_x_steps1]));
            }
        }
    }
}

#ifndef NDEBUG
template <typename POINT>
void Grid<POINT>::createGridGeometry(GeoLib::GEOObjects* geo_obj) const
{
    std::vector<std::string> grid_names;

    GeoLib::Point const& llf (getMinPoint());
    GeoLib::Point const& urb (getMaxPoint());

    const double dx ((urb[0] - llf[0]) / _n_steps[0]);
    const double dy ((urb[1] - llf[1]) / _n_steps[1]);
    const double dz ((urb[2] - llf[2]) / _n_steps[2]);

    // create grid names and grid boxes as geometry
    for (std::size_t i(0); i<_n_steps[0]; i++) {
        for (std::size_t j(0); j<_n_steps[1]; j++) {
            for (std::size_t k(0); k<_n_steps[2]; k++) {

                std::string name("Grid-");
                name += std::to_string(i);
                name += "-";
                name += std::to_string(j);
                name += "-";
                name += std::to_string (k);
                grid_names.push_back(name);

                {
                    auto points = std::unique_ptr<std::vector<GeoLib::Point*>>(
                        new std::vector<GeoLib::Point*>);
                    points->push_back(new GeoLib::Point(llf[0]+i*dx, llf[1]+j*dy, llf[2]+k*dz));
                    points->push_back(new GeoLib::Point(llf[0]+i*dx, llf[1]+(j+1)*dy, llf[2]+k*dz));
                    points->push_back(new GeoLib::Point(llf[0]+(i+1)*dx, llf[1]+(j+1)*dy, llf[2]+k*dz));
                    points->push_back(new GeoLib::Point(llf[0]+(i+1)*dx, llf[1]+j*dy, llf[2]+k*dz));
                    points->push_back(new GeoLib::Point(llf[0]+i*dx, llf[1]+j*dy, llf[2]+(k+1)*dz));
                    points->push_back(new GeoLib::Point(llf[0]+i*dx, llf[1]+(j+1)*dy, llf[2]+(k+1)*dz));
                    points->push_back(new GeoLib::Point(llf[0]+(i+1)*dx, llf[1]+(j+1)*dy, llf[2]+(k+1)*dz));
                    points->push_back(new GeoLib::Point(llf[0]+(i+1)*dx, llf[1]+j*dy, llf[2]+(k+1)*dz));
                    geo_obj->addPointVec(std::move(points), grid_names.back(),
                                         nullptr);
                }

                auto plys = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
                    new std::vector<GeoLib::Polyline*>);
                auto const& points = *geo_obj->getPointVec(grid_names.back());
                GeoLib::Polyline* ply0 (new GeoLib::Polyline(points));

                for (std::size_t l(0); l < 4; l++)
                    ply0->addPoint(l);
                ply0->addPoint(0);
                plys->push_back(ply0);

                GeoLib::Polyline* ply1 (new GeoLib::Polyline(points));
                for (std::size_t l(4); l < 8; l++)
                    ply1->addPoint(l);
                ply1->addPoint(4);
                plys->push_back(ply1);

                GeoLib::Polyline* ply2 (new GeoLib::Polyline(points));
                ply2->addPoint(0);
                ply2->addPoint(4);
                plys->push_back(ply2);

                GeoLib::Polyline* ply3 (new GeoLib::Polyline(points));
                ply3->addPoint(1);
                ply3->addPoint(5);
                plys->push_back(ply3);

                GeoLib::Polyline* ply4 (new GeoLib::Polyline(points));
                ply4->addPoint(2);
                ply4->addPoint(6);
                plys->push_back(ply4);

                GeoLib::Polyline* ply5 (new GeoLib::Polyline(points));
                ply5->addPoint(3);
                ply5->addPoint(7);
                plys->push_back(ply5);

                geo_obj->addPolylineVec(std::move(plys), grid_names.back(),
                    nullptr);
            }
        }
    }
    std::string merged_geo_name("Grid");

    geo_obj->mergeGeometries(grid_names, merged_geo_name);
}
#endif

template <typename POINT>
template <typename T>
std::array<std::size_t,3> Grid<POINT>::getGridCoords(T const& pnt) const
{
    std::array<std::size_t,3> coords;
    for (std::size_t k(0); k<3; k++) {
        if (pnt[k] < _min_pnt[k]) {
            coords[k] = 0;
        } else {
            if (pnt[k] > _max_pnt[k]) {
                coords[k] = _n_steps[k]-1;
            } else {
                coords[k] = static_cast<std::size_t>(
                    std::floor((pnt[k] - _min_pnt[k])) *
                    _inverse_step_sizes[k]);
            }
        }
    }
    return coords;
}

template <typename POINT>
template <typename P>
std::array<double,6> Grid<POINT>::getPointCellBorderDistances(P const& p,
    std::array<std::size_t,3> const& coords) const
{
    std::array<double,6> dists;
    dists[0] = std::abs(p[2]-_min_pnt[2] + coords[2]*_step_sizes[2]); // bottom
    dists[5] = std::abs(p[2]-_min_pnt[2] + (coords[2]+1)*_step_sizes[2]); // top

    dists[1] = std::abs(p[1]-_min_pnt[1] + coords[1]*_step_sizes[1]); // front
    dists[3] = std::abs(p[1]-_min_pnt[1] + (coords[1]+1)*_step_sizes[1]); // back

    dists[4] = std::abs(p[0]-_min_pnt[0] + coords[0]*_step_sizes[0]); // left
    dists[2] = std::abs(p[0]-_min_pnt[0] + (coords[0]+1)*_step_sizes[0]); // right
    return dists;
}

template <typename POINT>
template <typename P>
POINT* Grid<POINT>::getNearestPoint(P const& pnt) const
{
    std::array<std::size_t,3> coords(getGridCoords(pnt));

    double sqr_min_dist(MathLib::sqrDist(_min_pnt, _max_pnt));
    POINT* nearest_pnt(nullptr);

    std::array<double,6> dists(getPointCellBorderDistances(pnt, coords));

    if (calcNearestPointInGridCell(pnt, coords, sqr_min_dist, nearest_pnt)) {
        double min_dist(sqrt(sqr_min_dist));
        if (dists[0] >= min_dist && dists[1] >= min_dist
            && dists[2] >= min_dist && dists[3] >= min_dist
            && dists[4] >= min_dist && dists[5] >= min_dist) {
            return nearest_pnt;
        }
    } else {
        // search in all border cells for at least one neighbor
        double sqr_min_dist_tmp;
        POINT * nearest_pnt_tmp(nullptr);
        std::size_t offset(1);

        while (nearest_pnt == nullptr) {
            std::array<std::size_t,3> ijk{{
                coords[0]<offset ? 0 : coords[0]-offset,
                coords[1]<offset ? 0 : coords[1]-offset,
                coords[2]<offset ? 0 : coords[2]-offset}};
            for (; ijk[0]<coords[0]+offset; ijk[0]++) {
                for (; ijk[1] < coords[1] + offset; ijk[1]++) {
                    for (; ijk[2] < coords[2] + offset; ijk[2]++) {
                        // do not check the origin grid cell twice
                        if (ijk[0] == coords[0]
                            && ijk[1] == coords[1]
                            && ijk[2] == coords[2]) {
                            continue;
                        }
                        // check if temporary grid cell coordinates are valid
                        if (ijk[0] >= _n_steps[0]
                            || ijk[1] >= _n_steps[1]
                            || ijk[2] >= _n_steps[2]) {
                            continue;
                        }

                        if (calcNearestPointInGridCell(pnt, ijk,
                                sqr_min_dist_tmp, nearest_pnt_tmp)) {
                            if (sqr_min_dist_tmp < sqr_min_dist) {
                                sqr_min_dist = sqr_min_dist_tmp;
                                nearest_pnt = nearest_pnt_tmp;
                            }
                        }
                    } // end k
                } // end j
            } // end i
            offset++;
        } // end while
    } // end else

    double len(sqrt(MathLib::sqrDist(pnt, *nearest_pnt)));
    // search all other grid cells within the cube with the edge nodes
    std::vector<std::vector<POINT*> const*> vecs_of_pnts(
        getPntVecsOfGridCellsIntersectingCube(pnt, len));

    const std::size_t n_vecs(vecs_of_pnts.size());
    for (std::size_t j(0); j<n_vecs; j++) {
        std::vector<POINT*> const& pnts(*(vecs_of_pnts[j]));
        const std::size_t n_pnts(pnts.size());
        for (std::size_t k(0); k<n_pnts; k++) {
            const double sqr_dist(MathLib::sqrDist(pnt, *pnts[k]));
            if (sqr_dist < sqr_min_dist) {
                sqr_min_dist = sqr_dist;
                nearest_pnt = pnts[k];
            }
        }
    }

    return nearest_pnt;
}

template <typename POINT>
void Grid<POINT>::initNumberOfSteps(std::size_t n_per_cell,
    std::size_t n_pnts, std::array<double,3> const& extensions)
{
    double const max_extension(
        *std::max_element(extensions.cbegin(), extensions.cend()));

    std::bitset<3> dim; // all bits set to zero
    for (std::size_t k(0); k<3; ++k) {
        // set dimension if the ratio kth-extension/max_extension >= 1e-4
        if (extensions[k] >= 1e-4 * max_extension) {
            dim[k] = true;
        }
    }

    // structured grid: n_cells = _n_steps[0] * _n_steps[1] * _n_steps[2]
    // *** condition: n_pnts / n_cells < n_per_cell
    // => n_pnts / n_per_cell < n_cells
    // _n_steps[1] = _n_steps[0] * extensions[1]/extensions[0],
    // _n_steps[2] = _n_steps[0] * extensions[2]/extensions[0],
    // => n_cells = _n_steps[0]^3 * extensions[1]/extensions[0] *
    //              extensions[2]/extensions[0],
    // => _n_steps[0] = cbrt(n_cells * extensions[0]^2 /
    //                       (extensions[1]*extensions[2]))
    auto sc_ceil = [](double v) {
        return static_cast<std::size_t>(ceil(v));
    };

    switch (dim.count()) {
    case 3: // 3d case
        _n_steps[0] =
            sc_ceil(std::cbrt(n_pnts * (extensions[0] / extensions[1]) *
                              (extensions[0] / extensions[2]) / n_per_cell));
        _n_steps[1] = sc_ceil(_n_steps[0] *
                              std::min(extensions[1] / extensions[0], 100.0));
        _n_steps[2] = sc_ceil(_n_steps[0] *
                              std::min(extensions[2] / extensions[0], 100.0));
        break;
    case 2: // 2d cases
        if (dim[0] && dim[1]) { // xy
            _n_steps[0] = sc_ceil(std::sqrt(n_pnts * extensions[0] /
                                            (n_per_cell * extensions[1])));
            _n_steps[1] = sc_ceil(
                _n_steps[0] * std::min(extensions[1] / extensions[0], 100.0));
        } else if (dim[0] && dim[2]) { // xz
            _n_steps[0] = sc_ceil(std::sqrt(n_pnts * extensions[0] /
                                            (n_per_cell * extensions[2])));
            _n_steps[2] = sc_ceil(
                _n_steps[0] * std::min(extensions[2] / extensions[0], 100.0));
        } else if (dim[1] && dim[2]) { // yz
            _n_steps[1] = sc_ceil(std::sqrt(n_pnts * extensions[1] /
                                            (n_per_cell * extensions[2])));
            _n_steps[2] = sc_ceil(std::min(extensions[2]/extensions[1],100.0));
        }
        break;
    case 1: // 1d cases
        for (std::size_t k(0); k<3; ++k) {
            if (dim[k]) {
                _n_steps[k] = sc_ceil(static_cast<double>(n_pnts)/n_per_cell);
            }
        }
    }
}

template <typename POINT>
template <typename P>
bool Grid<POINT>::calcNearestPointInGridCell(P const& pnt,
    std::array<std::size_t,3> const& coords,
    double &sqr_min_dist,
    POINT* &nearest_pnt) const
{
    const std::size_t grid_idx (coords[0]+coords[1]*_n_steps[0]+coords[2]*_n_steps[0]*_n_steps[1]);
    std::vector<typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type> const& pnts(_grid_cell_nodes_map[grid_idx]);
    if (pnts.empty()) return false;

    const std::size_t n_pnts(pnts.size());
    sqr_min_dist = MathLib::sqrDist(*pnts[0], pnt);
    nearest_pnt = pnts[0];
    for (std::size_t i(1); i < n_pnts; i++) {
        const double sqr_dist(MathLib::sqrDist(*pnts[i], pnt));
        if (sqr_dist < sqr_min_dist) {
            sqr_min_dist = sqr_dist;
            nearest_pnt = pnts[i];
        }
    }
    return true;
}

template <typename POINT>
template <typename P>
std::vector<std::size_t>
Grid<POINT>::getPointsInEpsilonEnvironment(P const& pnt, double eps) const
{
    std::vector<std::vector<POINT*> const*> vec_pnts(
        getPntVecsOfGridCellsIntersectingCube(pnt, eps));

    double const sqr_eps(eps*eps);

    std::vector<std::size_t> pnts;
    for (auto vec : vec_pnts) {
        for (auto const p : *vec) {
            if (MathLib::sqrDist(*p, pnt) < sqr_eps) {
                pnts.push_back(p->getID());
            }
        }
    }

    return pnts;
}

} // end namespace GeoLib

#endif /* MESHGRID_H_ */
