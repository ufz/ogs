/**
 * \file
 * \author Thomas Fischer
 * \date   no date
 * \brief  Definition of sparse matrix IO functions.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include <cassert>

//extern void CS_write(char*, unsigned, unsigned const*, unsigned const*, double const*);
//extern void CS_read(char*, unsigned&, unsigned*&, unsigned*&, double*&);

template<class T> void CS_write(std::ostream &os, unsigned n, unsigned const* iA, unsigned const* jA, T const* A)
{
	os.write(reinterpret_cast<char*>(&n), sizeof(unsigned));
	os.write(reinterpret_cast<char*>(const_cast<unsigned*>(iA)), (n + 1) * sizeof(unsigned));
	os.write(reinterpret_cast<char*>(const_cast<unsigned*>(jA)), iA[n] * sizeof(unsigned));
	os.write(reinterpret_cast<char*>(A), iA[n] * sizeof(T));
}

template<class T> void CS_read(std::istream &is, unsigned &n, unsigned* &iA, unsigned* &jA, T* &A)
{
	is.read(reinterpret_cast<char*>(&n), sizeof(unsigned));
	if (iA != NULL) {
		delete[] iA;
		delete[] jA;
		delete[] A;
	}
	iA = new unsigned[n + 1];
	assert(iA != NULL);
	is.read(reinterpret_cast<char*>(iA), (n + 1) * sizeof(unsigned));

	jA = new unsigned[iA[n]];
	assert(jA != NULL);
	is.read(reinterpret_cast<char*>(jA), iA[n] * sizeof(unsigned));

	A = new T[iA[n]];
	assert(A != NULL);
	is.read(reinterpret_cast<char*>(A), iA[n] * sizeof(T));

#ifndef NDEBUG
	// do simple checks
	if (iA[0] != 0) std::cerr << "\nCRS matrix: array iA doesn't start with 0\n";

	unsigned i = 0;
	while (i < iA[n] && jA[i] < n)
		++i;
	if (i < iA[n]) std::cerr << "\nCRS matrix: the " << i
					<< "th entry of jA has the value " << jA[i] << ", which is out of bounds.\n";
#endif
}

#endif

