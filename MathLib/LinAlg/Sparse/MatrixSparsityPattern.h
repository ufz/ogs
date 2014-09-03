/**
 * @file MatrixSparsityPattern.h
 * @author Thomas Fischer
 * @date Apr 12, 2013
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef MATRIXSPARSITYPATTERN_H_
#define MATRIXSPARSITYPATTERN_H_

#include <set>
#include <vector>

#include "BaseLib/CodingTools.h"

namespace MathLib
{

/// \brief Class for representation of matrix sparsity pattern, required for
/// creation of sparse matrices.
///
/// \details Current implementation requires only number of matrix's rows,
/// allowing thus non-rectangular sparsity patterns. The class is based on
/// std::set container which automatically check for multiple insertions.
class MatrixSparsityPattern
{
public:
	/// Constant iterator over sorted entries of a row.
	typedef std::set<std::size_t>::const_iterator ConstRowIterator;

	explicit MatrixSparsityPattern(std::size_t const n_rows);
	virtual ~MatrixSparsityPattern();

	/// Returns number of sparsity pattern rows.
	std::size_t getNRows() const;

	/// Constant iterator over sorted entries of a row.
	ConstRowIterator getRowBeginIterator(std::size_t const row) const;
	/// Constant iterator over sorted entries of a row.
	ConstRowIterator getRowEndIterator(std::size_t const row) const;

	/// Inserts an entry in the sparsity pattern.
	/// \param row The row index the entry should be inserted to. The row index must be less or equal to the value returned by getNRows().
	/// \param col The column index. A new entry will be created if needed.
	void insert(std::size_t const row, std::size_t const col);

private:
	DISALLOW_COPY_AND_ASSIGN(MatrixSparsityPattern);

	std::vector<std::set<std::size_t> > _pattern;
};

} // end namespace MathLib

#endif // MATRIXSPARSITYPATTERN_H_
