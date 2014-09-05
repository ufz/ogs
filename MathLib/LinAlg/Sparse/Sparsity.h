/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-06-25
 * \brief  Definition of the Sparsity class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SPARSITY_H_
#define SPARSITY_H_

#include <vector>
#include <set>

namespace MathLib
{

/**
 * \brief Row-major sparse pattern
 */
typedef std::vector<std::set<std::size_t> > RowMajorSparsity;

} //MathLib

#endif // SPARSITY_H_
