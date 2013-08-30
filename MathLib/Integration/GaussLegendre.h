/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

#include <utility>

namespace MathLib
{

/**
 * \brief Gauss-Legendre quadrature method
 *
 */
class GaussLegendre
{
public:
    /**
     * return a point and weight
     *
     * \param n_points   the number of quadrature points
     * \param point_id   point index
     * \return point, weight at the given index
     */
    static std::pair<double, double> getPoint(unsigned n_points, unsigned point_id);

    /**
     * integrate a given function over [-1, 1]
     *
     * \tparam T    Function type
     * \param f     Function object
     * \param nQ    The number of sampling points
     */
    template <class T>
    static double integrate(T &f, unsigned nQ)
    {
        double val = .0;
        for (unsigned i=0; i<nQ; i++) {
            auto const& pt = getPoint(nQ, i);
            val += f(pt.first) * pt.second;
        }
        return val;
    }
};

}

#endif // GAUSSLEGENDRE_H_

