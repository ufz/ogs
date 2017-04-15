/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-24
 * \brief  Implementation of the StrictDoubleValidator class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <QDoubleValidator>

/**
 * \brief A validator for an input field which only accepts decimals.
 * Source code adapted from [Qt developer fac](http://developer.qt.nokia.com/faq/answer/i_can_still_insert_numbers_outside_the_range_specified_with_a_qdoublevalida)
 */
class StrictDoubleValidator : public QDoubleValidator
{
public:
    StrictDoubleValidator ( double min, double max, std::size_t decimals, QObject* parent = nullptr) :
        QDoubleValidator( min, max, decimals, parent)
    {}

    StrictDoubleValidator ( QObject* parent = nullptr) :
        QDoubleValidator( parent)
    {}

    QValidator::State validate(QString& input, int& pos) const override
    {
        if (input.isEmpty() || input == "." || input == "-") return Intermediate;

        if (QDoubleValidator::validate(input, pos) != Acceptable)
            return Invalid;
        return Acceptable;
    }
};
