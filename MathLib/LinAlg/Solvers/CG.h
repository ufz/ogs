/**
 * \file
 * \author Thomas Fischer
 * \date   2011-09-27
 * \brief  Definition of the CG functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CG_H_
#define CG_H_

namespace MathLib {

// forward declaration
template <typename PF_TYPE, typename IDX_TYPE> class CRSMatrix;

unsigned CG(CRSMatrix<double,unsigned> const * mat, double const * const b,
		double* const x, double& eps, unsigned& nsteps);

#ifdef _OPENMP
unsigned CGParallel(CRSMatrix<double,unsigned> const * mat, double const * const b,
		double* const x, double& eps, unsigned& nsteps);
#endif

} // end namespace MathLib

#endif /* SOLVER_H_ */
