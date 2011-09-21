#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include <cassert>

extern void CS_write(char*, unsigned, unsigned*, unsigned*, double*);
extern void CS_read(char*, unsigned&, unsigned*&, unsigned*&, double*&);

template<class T> void
CS_write(std::ostream &os, unsigned n, unsigned* iA, unsigned* jA, T* A)
{
  os.write((char*) &n, sizeof(unsigned));
  os.write((char*) iA, (n+1)*sizeof(unsigned));
  os.write((char*) jA, iA[n]*sizeof(unsigned));
  os.write((char*) A, iA[n]*sizeof(T));
}

template<class T> void
CS_read(std::istream &is, unsigned &n, unsigned* &iA, unsigned* &jA, T* &A)
{
  is.read((char*) &n, sizeof(unsigned));
  if (iA != NULL) {
    delete [] iA;
    delete [] jA;
    delete [] A;
  }
  iA = new unsigned[n+1];
  assert(iA!=NULL);
  is.read((char*) iA, (n+1)*sizeof(unsigned));

  jA = new unsigned[iA[n]];
  assert(jA!=NULL);
  is.read((char*) jA, iA[n]*sizeof(unsigned));

  A = new T[iA[n]];
  assert(A!=NULL);
  is.read((char*) A, iA[n]*sizeof(T));

#ifndef NDEBUG
  // do simple checks
  if (iA[0]!=0)
    std::cerr << std::endl
              << "CRS matrix: array iA doesn't start with 0" << std::endl;

  unsigned i = 0;
  while (i<iA[n] && jA[i]<n) ++i;
  if (i<iA[n])
    std::cerr << std::endl
              << "CRS matrix: the " << i << "th entry of jA has the value "
              << jA[i] << ", which is out of bounds." << std::endl;
#endif
}

void CS_transp(unsigned n, double* A, unsigned* jA, unsigned *iA,
               double* B, unsigned *jB, unsigned *iB)
{
  unsigned nnz = iA[n];
  unsigned *inz = new unsigned[n];

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

#endif
