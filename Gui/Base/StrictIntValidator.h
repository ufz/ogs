/*
 * StrictIntValidator.h
 *
 *  Created on: Jan 24, 2011
 *      Author: TF
 */

#ifndef STRICTINTVALIDATOR_H_
#define STRICTINTVALIDATOR_H_

#include <QIntValidator>

/**
 * \brief A validator for an input field which only accepts integers.
 * Source code adapted from [Qt developer fac](http://developer.qt.nokia.com/faq/answer/i_can_still_insert_numbers_outside_the_range_specified_with_a_qdoublevalida)
 */
class StrictIntValidator : public QIntValidator
{
public:
	StrictIntValidator ( int min, int max, QObject* parent = 0) :
		QIntValidator( min, max, parent)
	{}

	QValidator::State validate(QString & input, int &pos) const
	{
		if (input.isEmpty()) return Intermediate;

		if (QIntValidator::validate(input, pos) != Acceptable)
			return Invalid;
		return Acceptable;
	}
};

#endif /* STRICTINTVALIDATOR_H_ */
