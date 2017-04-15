/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cmath>

namespace BaseLib
{

/**
 * Interface class for subdivision operators
 */
class ISubdivision
{
public:
    /// Returns a vector of subdivided points
    virtual std::vector<double> operator()() const = 0;

    virtual ~ISubdivision() = default;
};

/**
 * Uniform subdivision operator
 */
class UniformSubdivision : public ISubdivision
{
public:
    /**
     * Configuration
     * @param length          total length to be subdivided
     * @param n_subdivision   the number of subdivision
     */
    UniformSubdivision(double length, std::size_t n_subdivision)
    : _length(length), _n_subdivision(n_subdivision) {}

    /// Returns a vector of subdivided points
    std::vector<double> operator()() const override
    {
        std::vector<double> x;
        x.reserve(_n_subdivision+1);
        const double dL = _length/static_cast<double>(_n_subdivision);
        for (std::size_t i=0; i<_n_subdivision+1; i++)
            x.push_back(i*dL);
        return x;
    }

private:
    const double _length;
    const std::size_t _n_subdivision;
};

/**
 * Gradual subdivision operator with a constant _multiplier
 */
class GradualSubdivision : public ISubdivision
{
public:
    /**
     * Constructor
     * @param L           total length to be subdivided
     * @param dL0         initial cell length
     * @param max_dL      maximum cell length
     * @param multiplier  multiplier to cell length
     */
    GradualSubdivision(
            const double L,
            const double dL0,
            const double max_dL,
            const double multiplier)
    : _length(L), _dL0(dL0), _max_dL(max_dL), _multiplier(multiplier) {}

    /// Returns a vector of subdivided points
    std::vector<double> operator()() const override
    {
        std::vector<double> vec_x;

        double x = 0;
        unsigned i=0;
        do {
            vec_x.push_back(x);
            x += std::min(_max_dL, _dL0*std::pow(_multiplier, static_cast<double>(i)));
            i++;
        } while (x<_length);

        if (vec_x.back() < _length) {
            double last_dx = vec_x[vec_x.size()-1] - vec_x[vec_x.size()-2];
            if (_length-vec_x.back()<last_dx)
                vec_x[vec_x.size()-1] = _length;
            else
                vec_x.push_back(_length);
        }
        return vec_x;
    }

private:
    const double _length;
    const double _dL0;
    const double _max_dL;
    const double _multiplier;
};

} // BaseLib
