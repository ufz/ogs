/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-24
 * \brief  Definition of the StrictIntValidator class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license

 */

#pragma once

#include <QIntValidator>

/**
 * \brief A validator for an input field which only accepts integers.
 * Source code adapted from [Qt developer fac](http://developer.qt.nokia.com/faq/answer/i_can_still_insert_numbers_outside_the_range_specified_with_a_qdoublevalida)
 */
class StrictIntValidator : public QIntValidator
{
public:
    StrictIntValidator(int min, int max, QObject* parent = nullptr)
        : QIntValidator(min, max, parent)
    {}

    QValidator::State validate(QString& input, int& pos) const override
    {
        if (input.isEmpty()) return Intermediate;

        if (QIntValidator::validate(input, pos) != Acceptable)
            return Invalid;
        return Acceptable;
    }
};
