/**
 * @file MatrixSparsityPattern.cpp
 * @author Thomas Fischer
 * @date Apr 12, 2013
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "LinAlg/Sparse/MatrixSparsityPattern.h"

namespace MathLib
{

MatrixSparsityPattern::MatrixSparsityPattern(std::size_t const n_rows) :
	_pattern(n_rows)
{}

MatrixSparsityPattern::~MatrixSparsityPattern()
{}

std::size_t MatrixSparsityPattern::getNRows() const
{
	return _pattern.size();
}

MatrixSparsityPattern::ConstRowIterator MatrixSparsityPattern::getRowBeginIterator(
        std::size_t const row) const
{
	return _pattern[row].cbegin();
}
MatrixSparsityPattern::ConstRowIterator MatrixSparsityPattern::getRowEndIterator(
        std::size_t const row) const
{
	return _pattern[row].cend();
}

void MatrixSparsityPattern::insert(std::size_t const row, std::size_t const col)
{
	_pattern[row].insert(col);
}

} // end namespace MathLib
