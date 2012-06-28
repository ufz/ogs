/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
 * \file generateDiagPrecond.h
 *
 * Created on 2011-09-28 by Thomas Fischer
 */

#ifndef GENERATEDIAGPRECOND_H_
#define GENERATEDIAGPRECOND_H_

namespace MathLib {

/**
 * diagonal preconditioner \f$P_{ii} = a_{ii}^{-1}\f$ associated with \f$n \times n\f$ matrix \f$A\f$
 * @param n number of rows / columns
 * @param iA row pointer of compressed row storage format
 * @param jA column index of compressed row storage format
 * @param A data entries of compressed row storage format
 * @param diag inverse entries of the diagonal
 * @return true, if all diagonal entries are distinct from zero, else false
 */
bool generateDiagPrecond(unsigned n, unsigned const*const iA, unsigned const*const jA,
				double const*const A, double* diag);

/**
 * diagonal preconditioner \f$P_{ii} = \left(\sum_{j} |a_{ij}|\right)^{-1}\f$ associated with \f$n \times n\f$ matrix \f$A\f$
 * @param n number of rows / columns
 * @param iA row pointer of compressed row storage format
 * @param A data entries of compressed row storage format
 * @param diag inverse entries of the diagonal
 * @return true, if all row sums are distinct from zero, else false
 */
bool generateDiagPrecondRowSum(unsigned n, unsigned const*const iA, double const*const A,
				double* diag);

/**
 * diagonal preconditioner \f$P_{ii} = \left(\max_{j} a_{ij}\right)^{-1}\f$ associated with \f$n \times n\f$ matrix \f$A\f$
 * @param n number of rows / columns
 * @param iA row pointer of compressed row storage format
 * @param A data entries of compressed row storage format
 * @param diag inverse entries of the diagonal
 * @return true, if all row sums are distinct from zero, else false
 */
bool generateDiagPrecondRowMax(unsigned n, unsigned const*const iA, double const*const A,
				double* diag);

} // end namespace MathLib

#endif /* PRECONDITIONER_H_ */
