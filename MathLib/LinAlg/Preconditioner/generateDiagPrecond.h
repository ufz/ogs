/*
 * generateDiagPrecond.h
 *
 *  Created on: Sep 28, 2011
 *      Author: TF
 */

#ifndef GENERATEDIAGPRECOND_H_
#define GENERATEDIAGPRECOND_H_

namespace MathLib {

bool generateDiagPrecond (unsigned n, unsigned* iA, unsigned* jA, double* A, double* diag);

} // end namespace MathLib

#endif /* PRECONDITIONER_H_ */
