/**
 * \author Norihiro Watanabe
 * \date   2013-08-30
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_FEM_ENUMS_H_
#define NUMLIB_FEM_ENUMS_H_

namespace NumLib
{

/**
  * \brief Shape matrix type to be calculated
  *
  */
enum class ShapeMatrixType
{
    N,      //< calculates N
    DNDR,   //< calculates dNdr
    N_J,    //< calculates N, dNdr, J, and detJ
    DNDR_J, //< caluclates dNdr, J, and detJ
    DNDX,   //< calculates dNdr, J, detJ, invJ, and dNdx
    ALL     //< calculates all
};

} // NumLib

#endif /* NUMLIB_FEM_ENUMS_H_ */
