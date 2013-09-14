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

#ifndef INTEGRATIONGAUSSREGULAR_H_
#define INTEGRATIONGAUSSREGULAR_H_

#include <cmath>
#include <array>

#include "MathLib/Integration/GaussLegendre.h"
#include "MathLib/TemplateWeightedPoint.h"


namespace NumLib
{

/**
 * \brief Gauss quadrature rule for regular shape elements, i.e. line, quad and hex
 *
 * \tpram N_DIM     Spatial dimension
 */
template <std::size_t N_DIM>
class IntegrationGaussRegular
{
public:
    /**
     * Construct this object with the given sampling level
     *
     * @param n_sampl_level     the sampling level (default 2)
     */
    explicit IntegrationGaussRegular(std::size_t n_sampl_level = 2)
    : _n_sampl_level(n_sampl_level), _n_sampl_pt(0)
    {
        this->setSamplingLevel(n_sampl_level);
    }

    /**
     * Change the sampling level
     *
     * @param n_sampl_level     the sampling level
     */
    void setSamplingLevel(std::size_t n_sampl_level)
    {
        this->_n_sampl_pt = std::pow(n_sampl_level, N_DIM);
        this->_n_sampl_level = n_sampl_level;
    }

    /// return current sampling level
    std::size_t getSamplingLevel() const {return _n_sampl_level;}

    /// return the number of sampling points
    std::size_t getNPoints() const {return _n_sampl_pt;}

    /**
     * get coordinates of a integration point
     *
     * @param igp      The integration point index
     * @return a weighted point
     */
    MathLib::TemplateWeightedPoint<double,double,N_DIM>
    getWeightedPoint(std::size_t igp) const
    {
        return getWeightedPoint(getSamplingLevel(), igp);
    }

    /**
     * get position indexes in r-s-t axis
     *
     * @param nGauss    The number of integration points
     * @param igp       The integration point index
     * @return  a tuple of position indexes
     */
    static std::array<std::size_t, N_DIM> getPosition(std::size_t nGauss, std::size_t igp);

    /**
     * get coordinates of a integration point
     *
     * @param nGauss    the number of integration points
     * @param pt_id     the sampling point id
     * @param x         coordinates
     * @return weight
     */
    static MathLib::TemplateWeightedPoint<double,double,N_DIM>
    getWeightedPoint(std::size_t nGauss, std::size_t igp);

private:
    /// Computes weighted point using given integration method.
    ///
    /// \tparam Method  Integration method to use.
    /// \param  pos     Point indices computed by getPosition.
    template <typename Method>
    static
    MathLib::TemplateWeightedPoint<double, double, N_DIM>
    getWeightedPoint(std::array<std::size_t, N_DIM> const& pos);

private:
    std::size_t _n_sampl_level;
    std::size_t _n_sampl_pt;
};

} // NumLib

#include "IntegrationGaussRegular.tpp"

#endif //INTEGRATIONGAUSSREGULAR_H_
