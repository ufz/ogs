/**
 * @file IsMatrixAssembled.h
 * @author Norihiro Watanabe
 * @date Jun 12, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef ISMATRIXASSEMBLED_H_
#define ISMATRIXASSEMBLED_H_

template <typename MAT_T>
inline bool isMatrixAssembled(MAT_T &)
{ return true; }

#endif /* ISMATRIXASSEMBLED_H_ */
