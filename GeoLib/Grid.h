/**
 * \file
 * \author Thomas Fischer
 * \date   2012-02-02
 * \brief  Definition of the Grid class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GRID_H_
#define GRID_H_

#include <vector>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// GeoLib
#include "AABB.h"
#include "GEOObjects.h"

#ifndef NDEBUG
// BaseLib
#include "StringTools.h"
#endif

// MathLib
#include "MathTools.h"

namespace GeoLib
{
template <typename POINT>
class Grid : public GeoLib::AABB<POINT>
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
	 * @param pnts (input) the points that are managed with the Grid
	 * @param max_num_per_grid_cell (input) max number per grid cell in the average (default 512)
	 *
	 */
	template <typename InputIterator>
	Grid(InputIterator first, InputIterator last, std::size_t max_num_per_grid_cell = 512) :
		GeoLib::AABB<POINT>(first, last),
		_n_steps({{1,1,1}}),
		_step_sizes({{0.0,0.0,0.0}}),
		_inverse_step_sizes({{0.0,0.0,0.0}}),
		_grid_cell_nodes_map(nullptr)
	{
		std::size_t n_pnts(std::distance(first, last));

		double delta[3] = { 0.0, 0.0, 0.0 };
		for (std::size_t k(0); k < 3; k++) {
			// make the bounding box a little bit bigger,
			// such that the node with maximal coordinates fits into the grid
			this->_max_pnt[k] += std::abs(this->_max_pnt[k]) * 1e-6;
			if (fabs(this->_max_pnt[k]) < std::numeric_limits<double>::epsilon()) {
				this->_max_pnt[k] = (this->_max_pnt[k] - this->_min_pnt[k]) * (1.0 + 1e-6);
			}
			delta[k] = this->_max_pnt[k] - this->_min_pnt[k];
		}

		// *** condition: n_pnts / (_n_steps[0] * _n_steps[1] * _n_steps[2]) < max_num_per_grid_cell
		// *** with _n_steps[1] = _n_steps[0] * delta[1]/delta[0], _n_steps[2] = _n_steps[0] * delta[2]/delta[0]
		const double eps(std::numeric_limits<double>::epsilon());
		if (fabs(delta[0]) < eps) { // dx == 0
			if(fabs(delta[1]) < eps) { // dy == 0
				if(fabs(delta[2]) < eps) { // degenerated case, dx == 0, dy == 0, dz == 0
					WARN("Grid constructor: Bounding volume [%f,%f] x [%f,%f] x [%f,%f] too small.",
					this->_min_pnt[0], this->_max_pnt[0],
					this->_min_pnt[1], this->_max_pnt[1],
					this->_min_pnt[2], this->_max_pnt[2]
					);
				} else { // 1d case: dx == 0, dy == 0, dz != 0
					_n_steps[0] = 1;
					_n_steps[1] = 1;
					_n_steps[2] = static_cast<std::size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
				}
			} else { // dy != 0
				if(fabs(delta[2]) < eps) { // 1d case: dx == 0, dy != 0, dz == 0
					_n_steps[0] = 1;
					_n_steps[1] = static_cast<std::size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
					_n_steps[2] = 1;
				} else { // 2d case: dx == 0, dy != 0, dz != 0
					_n_steps[0] = 1;
					_n_steps[1] = static_cast<std::size_t> (ceil(sqrt(n_pnts * delta[1] / (max_num_per_grid_cell * delta[2]))));
					_n_steps[2] = static_cast<std::size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
				}
			}
		} else { // dx != 0
			if(fabs(delta[1]) < eps) { // dy == 0
				if(fabs(delta[2]) < eps) { // 1d case: dx != 0, dy == 0, dz == 0
					_n_steps[0] = static_cast<std::size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
					_n_steps[1] = 1;
					_n_steps[2] = 1;
				} else { // 2d case: dx != 0, dy == 0, dz != 0
					_n_steps[0] = static_cast<std::size_t> (ceil(sqrt(n_pnts * delta[0] / (max_num_per_grid_cell * delta[2]))));
					_n_steps[1] = 1;
					_n_steps[2] = static_cast<std::size_t> (ceil(_n_steps[0] * delta[2] / delta[0]));
				}
			} else { // dy != 0
				if(fabs(delta[2]) < eps) { // 2d case: dx != 0, dy != 0, dz == 0
					_n_steps[0] = static_cast<std::size_t> (ceil(sqrt(n_pnts * delta[0] / (max_num_per_grid_cell * delta[1]))));
					_n_steps[1] = static_cast<std::size_t> (ceil(_n_steps[0] * delta[1] / delta[0]));
					_n_steps[2] = 1;
				} else { // 3d case: dx != 0, dy != 0, dz != 0
					_n_steps[0] = static_cast<std::size_t> (ceil(pow(n_pnts * delta[0] * delta[0]
									/ (max_num_per_grid_cell * delta[1] * delta[2]), 1. / 3.)));
					_n_steps[1] = std::max(static_cast<std::size_t>(1), std::min(static_cast<std::size_t> (ceil(_n_steps[0] * delta[1] / delta[0])), static_cast<std::size_t>(100)));
					_n_steps[2] = std::max(static_cast<std::size_t>(1), std::min(static_cast<std::size_t> (ceil(_n_steps[0] * delta[2] / delta[0])), static_cast<std::size_t>(100)));
				}
			}
		}

		const std::size_t n_plane(_n_steps[0] * _n_steps[1]);
		_grid_cell_nodes_map = new std::vector<POINT*>[n_plane * _n_steps[2]];

		// some frequently used expressions to fill the grid vectors
		for (std::size_t k(0); k < 3; k++) {
			if (fabs(delta[k]) < std::numeric_limits<double>::epsilon()) {
				delta[k] = std::numeric_limits<double>::epsilon();
			}
			_step_sizes[k] = delta[k] / _n_steps[k];
			_inverse_step_sizes[k] = 1.0 / _step_sizes[k];
		}

		// fill the grid vectors
		InputIterator it(first);
		while (it != last) {
			double const* const pnt(copyOrAddress(*it)->getCoords());
			const std::size_t i(static_cast<std::size_t> ((pnt[0] - this->_min_pnt[0]) * _inverse_step_sizes[0]));
			const std::size_t j(static_cast<std::size_t> ((pnt[1] - this->_min_pnt[1]) * _inverse_step_sizes[1]));
			const std::size_t k(static_cast<std::size_t> ((pnt[2] - this->_min_pnt[2]) * _inverse_step_sizes[2]));

			if (i < _n_steps[0] && j < _n_steps[1] && k < _n_steps[2]) {
				_grid_cell_nodes_map[i + j * _n_steps[0] + k * n_plane].push_back(const_cast<POINT*>(copyOrAddress(*it)));
			} else {
				ERR("Grid constructor: error computing indices [%d, %d, %d], max indices [%d, %d, %d].", i, j, k, _n_steps[0], _n_steps[1], _n_steps[2]);
			}
			it++;
		}

#ifndef NDEBUG
		std::size_t pnts_cnt(0);
		for (std::size_t k(0); k < n_plane * _n_steps[2]; k++)
			pnts_cnt += _grid_cell_nodes_map[k].size();

		assert(n_pnts==pnts_cnt);
#endif
	}

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
	 * POINT object a NULL pointer is returned.
	 *
	 * @param pnt a field that holds the coordinates of the point
	 * @return a pointer to the point with the smallest distance within the grid cells that are
	 * outlined above or NULL
	 */
	POINT* getNearestPoint(POINT const& pnt) const
	{
		std::size_t coords[3];
		getGridCoords(pnt, coords);

		double sqr_min_dist (MathLib::sqrDist(this->_min_pnt, this->_max_pnt));
		POINT* nearest_pnt(NULL);

		double dists[6];
		getPointCellBorderDistances(pnt, dists, coords);

		if (calcNearestPointInGridCell(pnt, coords, sqr_min_dist, nearest_pnt)) {
			double min_dist(sqrt(sqr_min_dist));
			if (dists[0] >= min_dist
							&& dists[1] >= min_dist
							&& dists[2] >= min_dist
							&& dists[3] >= min_dist
							&& dists[4] >= min_dist
							&& dists[5] >= min_dist) {
				return nearest_pnt;
			}
		} else {
			// search in all border cells for at least one neighbor
			double sqr_min_dist_tmp;
			POINT * nearest_pnt_tmp(NULL);
			std::size_t offset(1);

			while (nearest_pnt == NULL) {
				std::size_t tmp_coords[3];
				if (coords[0] < offset) {
					tmp_coords[0] = 0;
				} else {
					tmp_coords[0] = coords[0] - offset;
				}
				for (; tmp_coords[0] < coords[0] + offset; tmp_coords[0]++) {
					if (coords[1] < offset) {
						tmp_coords[1] = 0;
					} else {
						tmp_coords[1] = coords[1] - offset;
					}
					for (; tmp_coords[1] < coords[1] + offset; tmp_coords[1]++) {
						if (coords[2] < offset) {
							tmp_coords[2] = 0;
						} else {
							tmp_coords[2] = coords[2] - offset;
						}
						for (; tmp_coords[2] < coords[2] + offset; tmp_coords[2]++) {
							// do not check the origin grid cell twice
							if (!(tmp_coords[0] == coords[0] && tmp_coords[1] == coords[1]
							                && tmp_coords[2] == coords[2])) {
								// check if temporary grid cell coordinates are valid
								if (tmp_coords[0] < _n_steps[0] && tmp_coords[1] < _n_steps[1]
								                && tmp_coords[2] < _n_steps[2]) {
									if (calcNearestPointInGridCell(pnt, tmp_coords,
									                               sqr_min_dist_tmp,
									                               nearest_pnt_tmp)) {
										if (sqr_min_dist_tmp < sqr_min_dist) {
											sqr_min_dist = sqr_min_dist_tmp;
											nearest_pnt = nearest_pnt_tmp;
										}
									}
								} // valid grid cell coordinates
							} // same element
						} // end k
					} // end j
				} // end i
				offset++;
			} // end while
		} // end else

		double len (sqrt(MathLib::sqrDist(pnt.getCoords(), nearest_pnt->getCoords())));
		// search all other grid cells within the cube with the edge nodes
		std::vector<std::vector<POINT*> const*> vecs_of_pnts;
		getPntVecsOfGridCellsIntersectingCube(pnt, len, vecs_of_pnts);

		const std::size_t n_vecs(vecs_of_pnts.size());
		for (std::size_t j(0); j<n_vecs; j++) {
			std::vector<POINT*> const& pnts(*(vecs_of_pnts[j]));
			const std::size_t n_pnts(pnts.size());
			for (std::size_t k(0); k<n_pnts; k++) {
				const double sqr_dist (MathLib::sqrDist(pnt.getCoords(), pnts[k]->getCoords()));
				if (sqr_dist < sqr_min_dist) {
					sqr_min_dist = sqr_dist;
					nearest_pnt = pnts[k];
				}
			}
		}

		return nearest_pnt;
	}

	/**
	 * Method fetches the vectors of all grid cells intersecting the axis aligned cuboid
	 * defined by two points. The first point with minimal coordinates in all directions.
	 * The second point with maximal coordinates in all directions.
	 *
	 * @param center (input) the center point of the axis aligned cube
	 * @param half_len (input) half of the edge length of the axis aligned cube
	 * @param pnts (output) vector of vectors of points within grid cells that intersects
	 * the axis aligned cube
	 */
	void getPntVecsOfGridCellsIntersectingCube(POINT const& center, double half_len, std::vector<std::vector<POINT*> const*>& pnts) const;

	void getPntVecsOfGridCellsIntersectingCuboid(POINT const& min_pnt, POINT const& max_pnt, std::vector<std::vector<POINT*> const*>& pnts) const;


#ifndef NDEBUG
	/**
	 * Method creates a geometry for every mesh grid box. Additionally it
	 * creates one geometry containing all the box geometries.
	 * @param geo_obj
	 */
	void createGridGeometry(GeoLib::GEOObjects* geo_obj) const;
#endif

private:
	/**
	 * Method calculates the grid cell coordinates for the given point pnt. If
	 * the point is located outside of the bounding box the coordinates of the
	 * grid cell on the border that is nearest to given point will be returned.
	 * @param pnt (input) the coordinates of the point
	 * @param coords (output) the coordinates of the grid cell
	 */
	inline void getGridCoords(POINT const& pnt, std::size_t* coords) const;

	/**
	 *
	 * point numbering of the grid cell is as follow
	 * @code
	 *	     7 -------- 6
	 *	    /:         /|
	 *	   / :        / |
	 *	  /  :       /  |
	 *	 /   :      /   |
	 *	4 -------- 5    |
	 *	|    3 ....|... 2
	 *	|   .      |   /
	 *	|  .       |  /
	 *	| .        | /
	 *	|.         |/
	 *	0 -------- 1
	 *	@endcode
	 *	the face numbering is as follow:
	 *	face	nodes
	 *	0		0,3,2,1 bottom
	 *	1		0,1,5,4 front
	 *	2		1,2,6,5 right
	 *	3		2,3,7,6 back
	 *	4		3,0,4,7 left
	 *	5		4,5,6,7 top
	 * @param pnt (input) coordinates of the point
	 * @param dists (output) squared distances of the point to the faces
	 * ordered in the same sequence as above described
	 * @param coords coordinates of the grid cell
	 */
	void getPointCellBorderDistances(POINT const& pnt,
	                                 double dists[6],
	                                 std::size_t const* const coords) const;

	bool calcNearestPointInGridCell(POINT const& pnt, std::size_t const* const coords,
	                                double &sqr_min_dist,
	                                POINT* &nearest_pnt) const
	{
		const std::size_t grid_idx (coords[0] + coords[1] * _n_steps[0] + coords[2] * _n_steps[0] * _n_steps[1]);
		std::vector<typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type> const& pnts(_grid_cell_nodes_map[grid_idx]);
		if (pnts.empty()) return false;

		const std::size_t n_pnts(pnts.size());
		sqr_min_dist = MathLib::sqrDist(pnts[0]->getCoords(), pnt.getCoords());
		nearest_pnt = pnts[0];
		for (std::size_t i(1); i < n_pnts; i++) {
			const double sqr_dist(MathLib::sqrDist(pnts[i]->getCoords(), pnt.getCoords()));
			if (sqr_dist < sqr_min_dist) {
				sqr_min_dist = sqr_dist;
				nearest_pnt = pnts[i];
			}
		}
		return true;
	}

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

template<typename POINT>
void Grid<POINT>::getPntVecsOfGridCellsIntersectingCube(POINT const& center,
                                                        double half_len,
                                                        std::vector<std::vector<POINT*> const*>& pnts) const
{
	double tmp_pnt[3] = { center[0] - half_len, center[1] - half_len, center[2] - half_len }; // min
	std::size_t min_coords[3];
	getGridCoords(tmp_pnt, min_coords);

	tmp_pnt[0] = center[0] + half_len;
	tmp_pnt[1] = center[1] + half_len;
	tmp_pnt[2] = center[2] + half_len;
	std::size_t max_coords[3];
	getGridCoords(tmp_pnt, max_coords);

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

template<typename POINT>
void Grid<POINT>::getPntVecsOfGridCellsIntersectingCuboid(POINT const& min_pnt,
                                                          POINT const& max_pnt,
                                                          std::vector<std::vector<POINT*> const*>& pnts) const
{
	std::size_t min_coords[3];
	getGridCoords(min_pnt, min_coords);

	std::size_t max_coords[3];
	getGridCoords(max_pnt, max_coords);

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

	GeoLib::Point const& llf (this->getMinPoint());
	GeoLib::Point const& urb (this->getMaxPoint());

	const double dx ((urb[0] - llf[0]) / _n_steps[0]);
	const double dy ((urb[1] - llf[1]) / _n_steps[1]);
	const double dz ((urb[2] - llf[2]) / _n_steps[2]);

	// create grid names and grid boxes as geometry
	for (std::size_t i(0); i<_n_steps[0]; i++) {
		for (std::size_t j(0); j<_n_steps[1]; j++) {
			for (std::size_t k(0); k<_n_steps[2]; k++) {

				std::string name("Grid-");
				name += BaseLib::number2str<std::size_t>(i);
				name += "-";
				name += BaseLib::number2str<std::size_t>(j);
				name += "-";
				name += BaseLib::number2str<std::size_t> (k);
				grid_names.push_back(name);

				std::vector<GeoLib::Point*>* points (new std::vector<GeoLib::Point*>);
				points->push_back(new GeoLib::Point(llf[0]+i*dx, llf[1]+j*dy, llf[2]+k*dz));
				points->push_back(new GeoLib::Point(llf[0]+i*dx, llf[1]+(j+1)*dy, llf[2]+k*dz));
				points->push_back(new GeoLib::Point(llf[0]+(i+1)*dx, llf[1]+(j+1)*dy, llf[2]+k*dz));
				points->push_back(new GeoLib::Point(llf[0]+(i+1)*dx, llf[1]+j*dy, llf[2]+k*dz));
				points->push_back(new GeoLib::Point(llf[0]+i*dx, llf[1]+j*dy, llf[2]+(k+1)*dz));
				points->push_back(new GeoLib::Point(llf[0]+i*dx, llf[1]+(j+1)*dy, llf[2]+(k+1)*dz));
				points->push_back(new GeoLib::Point(llf[0]+(i+1)*dx, llf[1]+(j+1)*dy, llf[2]+(k+1)*dz));
				points->push_back(new GeoLib::Point(llf[0]+(i+1)*dx, llf[1]+j*dy, llf[2]+(k+1)*dz));
				geo_obj->addPointVec(points, grid_names[grid_names.size()-1], NULL);

				std::vector<GeoLib::Polyline*>* plys (
				        new std::vector<GeoLib::Polyline*>);
				GeoLib::Polyline* ply0 (new GeoLib::Polyline(*points));
				for (std::size_t l(0); l < 4; l++)
					ply0->addPoint(l);
				ply0->addPoint(0);
				plys->push_back(ply0);

				GeoLib::Polyline* ply1 (new GeoLib::Polyline(*points));
				for (std::size_t l(4); l < 8; l++)
					ply1->addPoint(l);
				ply1->addPoint(4);
				plys->push_back(ply1);

				GeoLib::Polyline* ply2 (new GeoLib::Polyline(*points));
				ply2->addPoint(0);
				ply2->addPoint(4);
				plys->push_back(ply2);

				GeoLib::Polyline* ply3 (new GeoLib::Polyline(*points));
				ply3->addPoint(1);
				ply3->addPoint(5);
				plys->push_back(ply3);

				GeoLib::Polyline* ply4 (new GeoLib::Polyline(*points));
				ply4->addPoint(2);
				ply4->addPoint(6);
				plys->push_back(ply4);

				GeoLib::Polyline* ply5 (new GeoLib::Polyline(*points));
				ply5->addPoint(3);
				ply5->addPoint(7);
				plys->push_back(ply5);

				geo_obj->addPolylineVec(plys,
				                        grid_names[grid_names.size() - 1], NULL);
			}
		}
	}
	std::string merged_geo_name("Grid");

	geo_obj->mergeGeometries(grid_names, merged_geo_name);
}
#endif

template <typename POINT>
void Grid<POINT>::getGridCoords(POINT const& pnt, std::size_t* coords) const
{
	for (std::size_t k(0); k<3; k++) {
		if (pnt[k] < this->_min_pnt[k]) {
			coords[k] = 0;
		} else {
			if (pnt[k] > this->_max_pnt[k]) {
				coords[k] = _n_steps[k]-1;
			} else {
				coords[k] = static_cast<std::size_t>((pnt[k]-this->_min_pnt[k]) * _inverse_step_sizes[k]);
			}
		}
	}
}

template <typename POINT>
void Grid<POINT>::getPointCellBorderDistances(POINT const& pnt,
                                              double dists[6],
                                              std::size_t const* const coords) const
{
	dists[0] = fabs(pnt[2] - this->_min_pnt[2] + coords[2] * _step_sizes[2]); // bottom
	dists[5] = fabs(pnt[2] - this->_min_pnt[2] + (coords[2] + 1) * _step_sizes[2]); // top

	dists[1] = fabs(pnt[1] - this->_min_pnt[1] + coords[1] * _step_sizes[1]); // front
	dists[3] = fabs(pnt[1] - this->_min_pnt[1] + (coords[1] + 1) * _step_sizes[1]); // back

	dists[4] = fabs(pnt[0] - this->_min_pnt[0] + coords[0] * _step_sizes[0]); // left
	dists[2] = fabs(pnt[0] - this->_min_pnt[0] + (coords[0] + 1) * _step_sizes[0]); // right
}
} // end namespace GeoLib

#endif /* MESHGRID_H_ */
