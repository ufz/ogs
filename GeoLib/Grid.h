/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Grid.h
 *
 * Created on 2012-02-02 by Thomas Fischer
 */

#ifndef GRID_H_
#define GRID_H_

#include <type_traits>
#include <vector>

// GeoLib
#include "AxisAlignedBoundingBox.h"
#include "GEOObjects.h"

#ifndef NDEBUG
// BaseLib
#include "StringTools.h"
#endif

namespace GeoLib {

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
	 * @param pnts (input) the points that are managed with the Grid
	 * @param max_num_per_grid_cell (input) max number per grid cell in the average (default 512)
	 *
	 * @note the somewhat wired template syntax chooses this constructor for non-pointer types
	 * at compile time
	 */
	template <typename T, typename = typename std::enable_if<std::is_same<T, POINT>::value && !std::is_pointer<T>::value>::type>
	Grid(std::vector<T> const& pnts, size_t max_num_per_grid_cell = 512) :
		GeoLib::AABB(), _grid_cell_nodes_map(NULL)
	{
		// compute axis aligned bounding box
		const size_t n_pnts(pnts.size());
		for (size_t k(0); k < n_pnts; k++) {
			this->update(pnts[k].getCoords());
		}

		double delta[3] = { 0.0, 0.0, 0.0 };
		for (size_t k(0); k < 3; k++) {
			// make the bounding box a little bit bigger,
			// such that the node with maximal coordinates fits into the grid
			_max_pnt[k] *= (1.0 + 1e-6);
			if (fabs(_max_pnt[k]) < std::numeric_limits<double>::epsilon()) {
				_max_pnt[k] = (_max_pnt[k] - _min_pnt[k]) * (1.0 + 1e-6);
			}
			delta[k] = _max_pnt[k] - _min_pnt[k];
		}

		// *** condition: n_pnts / (_n_steps[0] * _n_steps[1] * _n_steps[2]) < max_num_per_grid_cell
		// *** with _n_steps[1] = _n_steps[0] * delta[1]/delta[0], _n_steps[2] = _n_steps[0] * delta[2]/delta[0]
		if (fabs(delta[1]) < std::numeric_limits<double>::epsilon() ||
						fabs(delta[2]) < std::numeric_limits<double>::epsilon()) {
			// 1d case y = z = 0
			if (fabs(delta[1]) < std::numeric_limits<double>::epsilon() &&
							fabs(delta[2]) < std::numeric_limits<double>::epsilon()) {
				_n_steps[0] = static_cast<size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
				_n_steps[1] = 1;
				_n_steps[2] = 1;
			} else {
				// 1d case x = z = 0
				if (fabs(delta[0]) < std::numeric_limits<double>::epsilon() &&
								fabs(delta[2]) < std::numeric_limits<double>::epsilon()) {
					_n_steps[0] = 1;
					_n_steps[1] = static_cast<size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
					_n_steps[2] = 1;
				} else {
					// 1d case x = y = 0
					if (fabs(delta[0]) < std::numeric_limits<double>::epsilon() && fabs(delta[1])
									< std::numeric_limits<double>::epsilon()) {
						_n_steps[0] = 1;
						_n_steps[1] = 1;
						_n_steps[2] = static_cast<size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
					} else {
						// 2d case
						if (fabs(delta[1]) < std::numeric_limits<double>::epsilon()) {
							_n_steps[0] = static_cast<size_t> (ceil(sqrt(n_pnts * delta[0] / (max_num_per_grid_cell * delta[2]))));
							_n_steps[1] = 1;
							_n_steps[2] = static_cast<size_t> (ceil(_n_steps[0] * delta[2] / delta[0]));
						} else {
							_n_steps[0] = static_cast<size_t> (ceil(sqrt(n_pnts * delta[0] / (max_num_per_grid_cell * delta[1]))));
							_n_steps[1] = static_cast<size_t> (ceil(_n_steps[0] * delta[1] / delta[0]));
							_n_steps[2] = 1;
						}
					}
				}
			}
		} else {
			// 3d case
			_n_steps[0] = static_cast<size_t> (ceil(pow(n_pnts * delta[0] * delta[0]
							/ (max_num_per_grid_cell * delta[1] * delta[2]), 1. / 3.)));
			_n_steps[1] = static_cast<size_t> (ceil(_n_steps[0] * delta[1] / delta[0]));
			_n_steps[2] = static_cast<size_t> (ceil(_n_steps[0] * delta[2] / delta[0]));
		}

		const size_t n_plane(_n_steps[0] * _n_steps[1]);
		_grid_cell_nodes_map = new std::vector<POINT*>[n_plane * _n_steps[2]];

		// some frequently used expressions to fill the grid vectors
		for (size_t k(0); k < 3; k++) {
			_step_sizes[k] = delta[k] / _n_steps[k];
			_inverse_step_sizes[k] = 1.0 / _step_sizes[k];
		}

		// fill the grid vectors
		for (size_t l(0); l < n_pnts; l++) {
			double const* const pnt(pnts[l].getCoords());
			const size_t i(static_cast<size_t> ((pnt[0] - _min_pnt[0]) * _inverse_step_sizes[0]));
			const size_t j(static_cast<size_t> ((pnt[1] - _min_pnt[1]) * _inverse_step_sizes[1]));
			const size_t k(static_cast<size_t> ((pnt[2] - _min_pnt[2]) * _inverse_step_sizes[2]));

			if (i >= _n_steps[0] || j >= _n_steps[1] || k >= _n_steps[2]) {
				std::cout << "error computing indices " << std::endl;
			}

			_grid_cell_nodes_map[i + j * _n_steps[0] + k * n_plane].push_back(&pnts[l]);
		}

#ifndef NDEBUG
		size_t pnts_cnt(0);
		for (size_t k(0); k < n_plane * _n_steps[2]; k++)
			pnts_cnt += _grid_cell_nodes_map[k].size();

		assert(n_pnts==pnts_cnt);
#endif
	}

	/**
	 * @brief The constructor of the grid object takes a vector of pointers to points or nodes. This is
	 * the main difference to the previous constructor. Furthermore the
	 * user can specify the *average* maximum number of points per grid cell.
	 *
	 * The number of grid cells are computed with the following formula
	 * \f$\frac{n_{points}}{n_{cells}} \le n_{max\_per\_cell}\f$
	 *
	 * In order to limit the memory wasting the maximum number of points per grid cell
	 * (in the average) should be a power of two (since std::vector objects resize itself
	 * with this step size).
	 *
	 * @param pnts (input) a vector holding pointers to the points that are managed with the Grid
	 * @param max_num_per_grid_cell (input) max number per grid cell in the average (default 512)
	 *
	 * @note the somewhat wired template syntax chooses this constructor for pointer types at compile time
	 */
	template <typename T, bool dummy = true, typename std::enable_if<std::is_same<T, POINT>::value && std::is_pointer<T>::value, int>::type = 0>
	Grid(std::vector<T> const& pnts, size_t max_num_per_grid_cell = 512) :
		GeoLib::AABB(), _grid_cell_nodes_map(NULL)
	{
		// compute axis aligned bounding box
		const size_t n_pnts(pnts.size());
		for (size_t k(0); k < n_pnts; k++) {
			this->update(pnts[k]->getCoords());
		}

		double delta[3] = { 0.0, 0.0, 0.0 };
		for (size_t k(0); k < 3; k++) {
			// make the bounding box a little bit bigger,
			// such that the node with maximal coordinates fits into the grid
			_max_pnt[k] *= (1.0 + 1e-6);
			if (fabs(_max_pnt[k]) < std::numeric_limits<double>::epsilon()) {
				_max_pnt[k] = (_max_pnt[k] - _min_pnt[k]) * (1.0 + 1e-6);
			}
			delta[k] = _max_pnt[k] - _min_pnt[k];
		}

		// *** condition: n_pnts / (_n_steps[0] * _n_steps[1] * _n_steps[2]) < max_num_per_grid_cell
		// *** with _n_steps[1] = _n_steps[0] * delta[1]/delta[0], _n_steps[2] = _n_steps[0] * delta[2]/delta[0]
		if (fabs(delta[1]) < std::numeric_limits<double>::epsilon() ||
						fabs(delta[2]) < std::numeric_limits<double>::epsilon()) {
			// 1d case y = z = 0
			if (fabs(delta[1]) < std::numeric_limits<double>::epsilon() &&
							fabs(delta[2]) < std::numeric_limits<double>::epsilon()) {
				_n_steps[0] = static_cast<size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
				_n_steps[1] = 1;
				_n_steps[2] = 1;
			} else {
				// 1d case x = z = 0
				if (fabs(delta[0]) < std::numeric_limits<double>::epsilon() &&
								fabs(delta[2]) < std::numeric_limits<double>::epsilon()) {
					_n_steps[0] = 1;
					_n_steps[1] = static_cast<size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
					_n_steps[2] = 1;
				} else {
					// 1d case x = y = 0
					if (fabs(delta[0]) < std::numeric_limits<double>::epsilon() && fabs(delta[1])
									< std::numeric_limits<double>::epsilon()) {
						_n_steps[0] = 1;
						_n_steps[1] = 1;
						_n_steps[2] = static_cast<size_t> (ceil(n_pnts / (double) max_num_per_grid_cell));
					} else {
						// 2d case
						if (fabs(delta[1]) < std::numeric_limits<double>::epsilon()) {
							_n_steps[0] = static_cast<size_t> (ceil(sqrt(n_pnts * delta[0] / (max_num_per_grid_cell * delta[2]))));
							_n_steps[1] = 1;
							_n_steps[2] = static_cast<size_t> (ceil(_n_steps[0] * delta[2] / delta[0]));
						} else {
							_n_steps[0] = static_cast<size_t> (ceil(sqrt(n_pnts * delta[0] / (max_num_per_grid_cell * delta[1]))));
							_n_steps[1] = static_cast<size_t> (ceil(_n_steps[0] * delta[1] / delta[0]));
							_n_steps[2] = 1;
						}
					}
				}
			}
		} else {
			// 3d case
			_n_steps[0] = static_cast<size_t> (ceil(pow(n_pnts * delta[0] * delta[0]
							/ (max_num_per_grid_cell * delta[1] * delta[2]), 1. / 3.)));
			_n_steps[1] = static_cast<size_t> (ceil(_n_steps[0] * delta[1] / delta[0]));
			_n_steps[2] = static_cast<size_t> (ceil(_n_steps[0] * delta[2] / delta[0]));
		}

		const size_t n_plane(_n_steps[0] * _n_steps[1]);
		_grid_cell_nodes_map = new std::vector<typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type>[n_plane * _n_steps[2]];

		// some frequently used expressions to fill the grid vectors
		for (size_t k(0); k < 3; k++) {
			_step_sizes[k] = delta[k] / _n_steps[k];
			_inverse_step_sizes[k] = 1.0 / _step_sizes[k];
		}

		// fill the grid vectors
		for (size_t l(0); l < n_pnts; l++) {
			double const* const pnt(pnts[l]->getCoords());
			const size_t i(static_cast<size_t> ((pnt[0] - _min_pnt[0]) * _inverse_step_sizes[0]));
			const size_t j(static_cast<size_t> ((pnt[1] - _min_pnt[1]) * _inverse_step_sizes[1]));
			const size_t k(static_cast<size_t> ((pnt[2] - _min_pnt[2]) * _inverse_step_sizes[2]));

			if (i >= _n_steps[0] || j >= _n_steps[1] || k >= _n_steps[2]) {
				std::cout << "error computing indices " << std::endl;
			}

			_grid_cell_nodes_map[i + j * _n_steps[0] + k * n_plane].push_back(pnts[l]);
		}

#ifndef NDEBUG
		size_t pnts_cnt(0);
		for (size_t k(0); k < n_plane * _n_steps[2]; k++)
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
	const typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type
	getNearestPoint(double const*const pnt) const
	{
		size_t coords[3];
		getGridCoords(pnt, coords);

		double sqr_min_dist (MathLib::sqrDist(&_min_pnt, &_max_pnt));
		typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type nearest_pnt(NULL);

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
					typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type nearest_pnt_tmp(NULL);
					size_t offset(1);

					while (nearest_pnt == NULL) {
						size_t tmp_coords[3];
						for (tmp_coords[0] = coords[0]-offset; tmp_coords[0]<coords[0]+offset; tmp_coords[0]++) {
						for (tmp_coords[1] = coords[1]-offset; tmp_coords[1]<coords[1]+offset; tmp_coords[1]++) {
							for (tmp_coords[2] = coords[2]-offset; tmp_coords[2]<coords[2]+offset; tmp_coords[2]++) {
								// do not check the origin grid cell twice
								if (!(tmp_coords[0] == coords[0] && tmp_coords[1] == coords[1] && tmp_coords[2] == coords[2])) {
									// check if temporary grid cell coordinates are valid
									if (tmp_coords[0] < _n_steps[0] && tmp_coords[1] < _n_steps[1] && tmp_coords[2] < _n_steps[2]) {
										if (calcNearestPointInGridCell(pnt, tmp_coords, sqr_min_dist_tmp, nearest_pnt_tmp)) {
											if (sqr_min_dist_tmp < sqr_min_dist) {
												sqr_min_dist = sqr_min_dist_tmp;
												nearest_pnt = nearest_pnt_tmp;
											}
										}
									} // valid grid cell coordinates
								} // same element
							}  // end k
						} // end j
					} // end i
					offset++;
				} // end while
			} // end else

			double len (sqrt(MathLib::sqrDist(pnt, nearest_pnt->getCoords())));
			// search all other grid cells within the cube with the edge nodes
			std::vector<std::vector<typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type> const*> vecs_of_pnts;
			getVecsOfGridCellsIntersectingCube(pnt, len, vecs_of_pnts);

			const size_t n_vecs(vecs_of_pnts.size());
			for (size_t j(0); j<n_vecs; j++) {
				std::vector<typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type> const& pnts(*(vecs_of_pnts[j]));
				const size_t n_pnts(pnts.size());
				for (size_t k(0); k<n_pnts; k++) {
					const double sqr_dist (MathLib::sqrDist(pnt, pnts[k]->getCoords()));
					if (sqr_dist < sqr_min_dist) {
						sqr_min_dist = sqr_dist;
						nearest_pnt = pnts[k];
					}
				}
			}

			return nearest_pnt;
	}

	/**
	 * Method fetches the vectors of all grid cells intersecting the axis aligned cube
	 * defined by its center and half edge length.
	 *
	 * @param pnt (input) the center point of the axis aligned cube
	 * @param half_len (input) half of the edge length of the axis aligned cube
	 * @param pnts (output) vector of vectors of points within grid cells that intersects
	 * the axis aligned cube
	 */
	void getVecsOfGridCellsIntersectingCube(double const*const pnt, double half_len, std::vector<std::vector<POINT> const*>& pnts) const;

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
	inline void getGridCoords(double const*const pnt, size_t* coords) const;

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
	void getPointCellBorderDistances(double const*const pnt, double dists[6], size_t const* const coords) const;

	bool calcNearestPointInGridCell(double const* const pnt, size_t const* const coords,
					double &sqr_min_dist,
					typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type &nearest_pnt) const
	{
		const size_t grid_idx (coords[0] + coords[1] * _n_steps[0] + coords[2] * _n_steps[0] * _n_steps[1]);
		std::vector<typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type> const& pnts(_grid_cell_nodes_map[grid_idx]);
		if (pnts.empty()) return false;

		const size_t n_pnts(pnts.size());
		sqr_min_dist = MathLib::sqrDist(pnts[0]->getCoords(), pnt);
		nearest_pnt = pnts[0];
		for (size_t i(1); i < n_pnts; i++) {
			const double sqr_dist(MathLib::sqrDist(pnts[i]->getCoords(), pnt));
			if (sqr_dist < sqr_min_dist) {
				sqr_min_dist = sqr_dist;
				nearest_pnt = pnts[i];
			}
		}
		return true;
	}
	double _step_sizes[3];
	double _inverse_step_sizes[3];
	size_t _n_steps[3];
	/**
	 * This is an array that stores pointers to POINT objects. If POINT is a pointer type,
	 * std::remove_pointer returns the "base" type, else the POINT type is returned. Then,
	 * the base type is converted to a pointer type emploing std::add_pointer.
	 */
	std::vector<typename std::add_pointer<typename std::remove_pointer<POINT>::type>::type>* _grid_cell_nodes_map;
};

template<typename POINT>
void Grid<POINT>::getVecsOfGridCellsIntersectingCube(double const* const pnt, double half_len,
				std::vector<std::vector<POINT> const*>& pnts) const
{
	double tmp_pnt[3] = { pnt[0] - half_len, pnt[1] - half_len, pnt[2] - half_len }; // min
	size_t min_coords[3];
	getGridCoords(tmp_pnt, min_coords);

	tmp_pnt[0] = pnt[0] + half_len;
	tmp_pnt[1] = pnt[1] + half_len;
	tmp_pnt[2] = pnt[2] + half_len;
	size_t max_coords[3];
	getGridCoords(tmp_pnt, max_coords);

	size_t coords[3], steps0_x_steps1(_n_steps[0] * _n_steps[1]);
	for (coords[0] = min_coords[0]; coords[0] < max_coords[0] + 1; coords[0]++) {
		for (coords[1] = min_coords[1]; coords[1] < max_coords[1] + 1; coords[1]++) {
			const size_t coords0_p_coords1_x_steps0(coords[0] + coords[1] * _n_steps[0]);
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

	const double dx ((urb[0]-llf[0])/_n_steps[0]);
	const double dy ((urb[1]-llf[1])/_n_steps[1]);
	const double dz ((urb[2]-llf[2])/_n_steps[2]);

	// create grid names and grid boxes as geometry
	for (size_t i(0); i<_n_steps[0]; i++) {
		for (size_t j(0); j<_n_steps[1]; j++) {
			for (size_t k(0); k<_n_steps[2]; k++) {

				std::string name("Grid-");
				name += number2str<size_t>(i);
				name +="-";
				name += number2str<size_t>(j);
				name += "-";
				name += number2str<size_t> (k);
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

				std::vector<GeoLib::Polyline*>* plys (new std::vector<GeoLib::Polyline*>);
				GeoLib::Polyline* ply0 (new GeoLib::Polyline(*points));
				for (size_t l(0); l<4; l++) {
					ply0->addPoint(l);
				}
				ply0->addPoint(0);
				plys->push_back(ply0);

				GeoLib::Polyline* ply1 (new GeoLib::Polyline(*points));
				for (size_t l(4); l<8; l++) {
					ply1->addPoint(l);
				}
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

				geo_obj->addPolylineVec(plys, grid_names[grid_names.size()-1], NULL);
			}
		}
	}
	std::string merged_geo_name("Grid");

	geo_obj->mergeGeometries(grid_names, merged_geo_name);
}
#endif

template <typename POINT>
void Grid<POINT>::getGridCoords(double const*const pnt, size_t* coords) const
{
	for (size_t k(0); k<3; k++) {
		if (pnt[k] < _min_pnt[k]) {
			coords[k] = 0;
		} else {
			if (pnt[k] > _max_pnt[k]) {
				coords[k] = _n_steps[k]-1;
			} else {
				coords[k] = static_cast<size_t>((pnt[k]-_min_pnt[k]) * _inverse_step_sizes[k]);
			}
		}
	}
}

template <typename POINT>
void Grid<POINT>::getPointCellBorderDistances(double const*const pnt, double dists[6], size_t const* const coords) const
{

	dists[0] = (pnt[2] - _min_pnt[2] + coords[2]*_step_sizes[2]); // bottom
	dists[5] = (_step_sizes[2] - dists[0]); // top

	dists[1] = (pnt[1] - _min_pnt[1] + coords[1]*_step_sizes[1]); // front
	dists[3] = (_step_sizes[1] - dists[1]); // back

	dists[4] = (pnt[0] - _min_pnt[0] + coords[0]*_step_sizes[0]); // left
	dists[2] = (_step_sizes[0] - dists[4]); // right
}

} // end namespace GeoLib

#endif /* MESHGRID_H_ */
