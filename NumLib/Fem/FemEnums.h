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

/**
  * Bit flags representing whether shape calculation cares about a specific field or not.  1 = care, 0 = don't care.
  */
typedef enum
{
    SHAPE_N     =   0x01,
//    SHAPE_DNDR  =   0x02,
    SHAPE_DNDX  =   0x04,
//    SHAPE_J     =   0x08,
//    SHAPE_INV_J =   0x10,
    SHAPE_ALL   =   0x7f
} ShapeFieldType;

#endif /* NUMLIB_FEM_ENUMS_H_ */
