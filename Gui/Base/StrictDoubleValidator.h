/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-24
 * \brief  Implementation of the StrictDoubleValidator class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef STRICTDOUBLEVALIDATOR_H_
#define STRICTDOUBLEVALIDATOR_H_

#include <QDoubleValidator>

/**
 * \brief A validator for an input field which only accepts decimals.
 * Source code adapted from [Qt developer fac](http://developer.qt.nokia.com/faq/answer/i_can_still_insert_numbers_outside_the_range_specified_with_a_qdoublevalida)
 */
class StrictDoubleValidator : public QDoubleValidator
{
public:
	StrictDoubleValidator ( double min, double max, std::size_t decimals, QObject* parent = 0) :
		QDoubleValidator( min, max, decimals, parent)
	{}

	QValidator::State validate(QString & input, int &pos) const
	{
		if (input.isEmpty() || input == ".") return Intermediate;

		if (QDoubleValidator::validate(input, pos) != Acceptable)
			return Invalid;
		return Acceptable;
	}
};

#endif /* STRICTDOUBLEVALIDATOR_H_ */
