/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef EXTRAPOLATIONGAUSSQUAD_H_
#define EXTRAPOLATIONGAUSSQUAD_H_

#include <cstddef>

namespace NumLib
{

/**
 * \brief Extrapolation of Gauss point values to nodal values
 */
class ExtrapolationGaussQuad
{
public:
    /**
     * extrapolate gauss point values to nodal values
     *
     * @param gp_values     a vector of gauss point values
     * @param nodal_values  a vector of extrapolated nodal values
     */
    template <class T_GP_VEC, class T_NOD_VEC>
    static void extrapolate(const T_GP_VEC &gp_values, T_NOD_VEC &nodal_values);

private:
    static std::size_t getNodeIndexOfGaussQuad(std::size_t nGaussLevel, std::size_t igp);
    static double calculateXi_p(std::size_t nGaussLevel);
    static void getExtrapolatedPoints(std::size_t nodeIdOfGaussQuad, double Xi_p, double* pt_in_natural);

};

} // end namespace

#include "ExtrapolationGaussQuad.tpp"

#endif //EXTRAPOLATIONGAUSSQUAD_H_
