/**
 * \file
 * \author Thomas Fischer
 * \date   2011-10-04
 * \brief  Definition of the BiCGStab function.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BICGSTAB_H_
#define BICGSTAB_H_

namespace MathLib
{
template <typename FP_TYPE, typename IDX_TYPE>
class CRSMatrix;
}

namespace MathLib {

unsigned BiCGStab(CRSMatrix<double, unsigned> const& A, double* const b, double* const x,
                  double& eps, unsigned& nsteps);

} // end namespace MathLib

#endif /* BICGSTAB_H_ */
