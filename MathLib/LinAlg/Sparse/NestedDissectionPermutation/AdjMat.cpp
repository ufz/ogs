/**
 * \file AdjMat.cpp
 *
 *  Created on 2012-01-02 by Thomas Fischer
 */

// Base
#include "swap.h"
#include "quicksort.h"

#include "LinAlg/Sparse/NestedDissectionPermutation/AdjMat.h"

namespace MathLib {

AdjMat* AdjMat::getMat(unsigned beg, unsigned end,
		const unsigned* const op_perm, const unsigned* const po_perm) const
{
	const unsigned nsize(end - beg); // size of new matrix
	unsigned i, c; // row and col idx in permuted matrix
	unsigned j, idx; // pointer in jA
	unsigned r; // row idx in original matrix

	unsigned *iAn(new unsigned[nsize + 1]);
	iAn[0] = 0;

	unsigned *pos(new unsigned[nsize + 1]);
	for (i = 0; i <= nsize; i++)
		pos[i] = 0;

	for (i = beg; i < end; i++) {
		r = op_perm[i];
		idx = _row_ptr[r + 1];
		for (j = _row_ptr[r]; j < idx; j++) {
			c = po_perm[_col_idx[j]];
			if (beg <= c && c < end)
				++pos[i - beg];
		}
	}
	for (i = 0; i < nsize; i++)
		iAn[i + 1] = iAn[i] + pos[i];
	for (i = 0; i < nsize; i++)
		pos[i] = iAn[i];

	unsigned *jAn(new unsigned[iAn[nsize]]);
	for (i = beg; i < end; i++) {
		r = op_perm[i];
		idx = _row_ptr[r + 1];
		for (j = _row_ptr[r]; j < idx; j++) {
			c = po_perm[_col_idx[j]];
			if (beg <= c && c < end)
				jAn[pos[i - beg]++] = c - beg;
		}
	}

	delete[] pos;
	for (i = 0; i < nsize; ++i)
		BaseLib::quickSort(jAn, iAn[i], iAn[i + 1]);
	return new AdjMat(nsize, iAn, jAn, NULL);
}

/**
 * generate an adjacency matrix (the upper triangle part)
 * @param n number of nodes of matrix graph / number of rows/columns of the adjacency matrix
 * @param iA array of size of the number of rows/columns + 1, array contains pointer into jA array
 * @param jA array of the length of the number of non-zero entries (edges in the matrix graph)
 */
void genAdjMat(unsigned n, unsigned* &iA, unsigned* &jA)
{
	unsigned i;
	// count entries of each row
	unsigned* iAn = new unsigned[n + 1];
	for (i = 0; i <= n; ++i)
		iAn[i] = 0;

	// go through all strictly lower triangular entries (i,j) and check
	// whether (j,i) exists in the upper triangular part

	// set n pointers to the beginning of each row
	unsigned* co = new unsigned[n];
	for (i = 0; i < n; ++i)
		co[i] = iA[i];

	for (i = 0; i < n; ++i)
		for (unsigned k = iA[i]; k < iA[i + 1]; ++k) {
			unsigned j = jA[k];
			if (i < j)
				++iAn[i + 1]; // upper triangular entries count
			else { // lower triangular only if there is no counter part
				unsigned k1 = iA[j], k2 = iA[j + 1];
				if (i < jA[k1] || i > jA[k2 - 1])
					++iAn[j + 1]; // i is out of bounds
				else { // go through all uninspected entries in the jth row
					while (co[j] < k2 && i > jA[co[j]])
						++co[j];
					if (co[j] == k2 || i < jA[co[j]])
						++iAn[j + 1];
				}
			}
		}

	// construct array iAn by summing up the contents of iAn
	// con is a set of pointer refering to iAn
	unsigned* con = new unsigned[n];
	co[0] = con[0] = 0;
	for (i = 1; i < n; ++i) {
		co[i] = iA[i];
		con[i] = iAn[i];
		iAn[i + 1] += iAn[i];
	}

	unsigned *jAn = new unsigned[iAn[n]];
	for (i = 1; i < n; ++i)
		for (unsigned k = iA[i]; k < iA[i + 1]; ++k) {
			unsigned j = jA[k];
			// copy all transposed lower triangular entries and all upper
			// triangular elements up to that position
			if (j < i) {
				while (co[j] < iA[j + 1] && i > jA[co[j]]) {
					if (jA[co[j]] > j)
						jAn[con[j]++] = jA[co[j]];
					++co[j];
				}

				if (co[j] == iA[j + 1] || i <= jA[co[j]]) {
					jAn[con[j]++] = i;
					++co[i];
					if (i == jA[co[j]])
						++co[j];
				}
			}
		}

	// finish rows
	for (i = 0; i < n; ++i)
		for (unsigned k = co[i]; k < iA[i + 1]; ++k)
			if (i < jA[k])
				jAn[con[i]++] = jA[k];

	BaseLib::swap(jA, jAn);
	BaseLib::swap(iA, iAn);

	delete[] jAn;
	delete[] con;
	delete[] co;
	delete[] iAn;
}

void genFullAdjMat(unsigned n, unsigned* &iA, unsigned* &jA)
{
	unsigned i;
	// count entries of each column
	unsigned* cnt = new unsigned[n];
	for (i = 0; i < n; ++i)
		cnt[i] = 0;

	for (i = 0; i < n; ++i) {
		unsigned j = iA[i], idx = iA[i + 1];
		while (j < idx) {
			cnt[jA[j]]++;
			j++;
		}
	}

	// summing up entries
	for (i = 2; i < n; ++i)
		cnt[i] += cnt[i - 1];

	unsigned* iAn = new unsigned[n + 1]; // VALGRIND meldet hier Fehler
	iAn[0] = 0;
	for (i = 1; i <= n; ++i)
		iAn[i] = iA[i] + cnt[i - 1];

	unsigned *jAn = new unsigned[iAn[n]];
	for (unsigned k = 0; k < n; k++)
		cnt[k] = iAn[k];

	for (i = 0; i < n; ++i) {
		unsigned j = iA[i], idx = iA[i + 1];
		while (j < idx) {
			jAn[cnt[i]++] = jA[j];
			jAn[cnt[jA[j]]++] = i;
			j++;
		}
	}

	BaseLib::swap(jA, jAn);
	BaseLib::swap(iA, iAn);

	delete[] jAn;
	delete[] iAn;
	delete[] cnt;
}

void AdjMat::makeSymmetric()
{
	// store upper triangular mat values
	genAdjMat(MatrixBase::_n_rows, _row_ptr, _col_idx);
	// mirror the upper triangular part into lower
	genFullAdjMat(MatrixBase::_n_rows, _row_ptr, _col_idx);
}

} // end namespace MathLib
