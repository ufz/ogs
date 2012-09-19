/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file sparse.h
 *
 * Created on xxxx-xx-xx by Thomas Fischer
 */

#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include <cassert>

//extern void CS_write(char*, unsigned, unsigned const*, unsigned const*, double const*);
//extern void CS_read(char*, unsigned&, unsigned*&, unsigned*&, double*&);

template<class T> void CS_write(std::ostream &os, unsigned n, unsigned const* iA, unsigned const* jA, T const* A)
{
	os.write((char*) &n, sizeof(unsigned));
	os.write((char*) iA, (n + 1) * sizeof(unsigned));
	os.write((char*) jA, iA[n] * sizeof(unsigned));
	os.write((char*) A, iA[n] * sizeof(T));
}

template<class T> void CS_read(std::istream &is, unsigned &n, unsigned* &iA, unsigned* &jA, T* &A)
{
	is.read((char*) &n, sizeof(unsigned));
	if (iA != NULL) {
		delete[] iA;
		delete[] jA;
		delete[] A;
	}
	iA = new unsigned[n + 1];
	assert(iA != NULL);
	is.read((char*) iA, (n + 1) * sizeof(unsigned));

	jA = new unsigned[iA[n]];
	assert(jA != NULL);
	is.read((char*) jA, iA[n] * sizeof(unsigned));

	A = new T[iA[n]];
	assert(A != NULL);
	is.read((char*) A, iA[n] * sizeof(T));

#ifndef NDEBUG
	// do simple checks
	if (iA[0] != 0) std::cerr << std::endl << "CRS matrix: array iA doesn't start with 0"
					<< std::endl;

	unsigned i = 0;
	while (i < iA[n] && jA[i] < n)
		++i;
	if (i < iA[n]) std::cerr << std::endl << "CRS matrix: the " << i
					<< "th entry of jA has the value " << jA[i] << ", which is out of bounds."
					<< std::endl;
#endif
}

#endif

