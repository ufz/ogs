/**
 * \file
 * \author Thomas Fischer
 * \date   no date
 * \brief  Definition of vector IO functions.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VECTOR_IO_H
#define VECTOR_IO_H

#include <iostream>
#include <fstream>
#include <string>

// von http://www.zfx.info/Display/Thread.php?TID=177779
#include <sstream>

#ifdef UNICODE
typedef std::wstringstream unistringstream;
typedef std::wstring unistring;
#else
typedef std::stringstream unistringstream;
typedef std::string unistring;
#endif

template<typename A, typename T> inline const A lexical_cast(const T& source)
{
	unistringstream s;

	s << source;

	A destination;
	s >> destination;

	return (destination);
}

/** reads the number of lines of non-binary stream */
unsigned readLines ( std::istream &in )
{
	unsigned k = 0;
	while (!in.eof()) {
		std::string str;
		getline (in, str);
		k++;
	}

	return k;
}

/** reads N values of non-binary stream */
template <class T> void read ( std::istream &in, unsigned N, T *a )
{
	unsigned k = 0;
	std::string ws (" \t");

	while (!in.eof () and k <= N) {
		std::string t;
		 getline (in, t);
		std::size_t i1;
		if (t.length () != 0) {
			i1 = t.find_first_not_of (ws, 0);
			if (i1 != std::string::npos) {
				std::size_t i2 = t.find_first_of (ws, i1);
				if (i2 != std::string::npos) {
					a[k++] = lexical_cast<T> ( t.substr(i1, i2-i1) );
				} else {
					a[k++] = lexical_cast<T> ( t.substr(i1, t.size()-i1));
				}
			}
		}
	}
}

template <class T> int readArrays ( const std::string &fname, std::size_t &s,
	std::size_t n, T* &arr )
{
	// open stream
	std::ifstream in (fname.c_str());
	// read number of rows
	if (in) {
		s = readLines(in)-1;
		in.close();
	} else {
		std::cout << "could not open " << fname << std::endl;
		return 1;
	}
	if (s > 0) {
		arr = new T[s * n];
		std::size_t j(0), l(0);
		// read values
		in.open (fname.c_str(), std::ios::in);
		for (j=0; j<s; ++j) {
			for (l=0; l<n; l++) in >> arr[l*s+j];
		}
		in.close ();
	}
	return 0;
}

#endif

