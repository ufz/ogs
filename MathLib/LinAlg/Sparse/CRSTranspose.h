/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file CRSTranspose.h
 *
 * Created on 2011-09-27 by Thomas Fischer
 */

#ifndef CRSTRANSPOSE_H_
#define CRSTRANSPOSE_H_

void CS_transp(unsigned n, unsigned *iA, unsigned* jA, double* A,
				unsigned *iB, unsigned *jB, double* B)
{
  unsigned nnz = iA[n];
  unsigned *inz(new unsigned[n]);

  for (unsigned i=0; i<n; ++i) inz[i] = 0;

  // compute number of entries of each column in A
  for (unsigned l=0; l<nnz; ++l) inz[jA[l]]++;

  // create iB
  iB[0] = nnz = 0;
  for (unsigned i=0; i<n; ++i) {
    nnz += inz[i];
    inz[i] = 0;
    iB[i+1] = nnz;
  }

  // create arrays jB, B
  unsigned l = iA[0];
  for (unsigned i=0; i<n; ++i) {
    while (l<iA[i+1]) {
      unsigned j = jA[l];
      unsigned k = iB[j] + inz[j]++;
      B[k] = A[l++];
      jB[k] = i;
    }
  }

  delete [] inz;
}

#endif /* CRSTRANSPOSE_H_ */
