/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Sparsity.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
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

} // end namespace MathLib

#endif // SPARSITY_H_
