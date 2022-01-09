/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-24
 * \brief  Implementation of the StrictDoubleValidator class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <QDoubleValidator>

/**
 * \brief A validator for an input field which only accepts decimals.
 * Source code adapted from [StackOverflow](https://stackoverflow.com/questions/19571033/allow-entry-in-qlineedit-only-within-range-of-qdoublevalidator)
 */
class StrictDoubleValidator : public QDoubleValidator
{
public:
    StrictDoubleValidator ( double min, double max, std::size_t decimals, QObject* parent = nullptr) :
        QDoubleValidator( min, max, decimals, parent)
    {}

    explicit StrictDoubleValidator(QObject* parent = nullptr)
        : QDoubleValidator(parent)
    {}

    QValidator::State validate(QString& input, int& pos) const override
    {
        Q_UNUSED(pos);
        if (input.isEmpty() || input == "." || input == "-") return Intermediate;

        QChar const decimalPoint('.');
        if (input.indexOf(decimalPoint) != -1)
        {
            int const charsAfterPoint = input.length() - input.indexOf(decimalPoint) - 1;
            if (charsAfterPoint > decimals())
                return QValidator::Invalid;
        }

        bool ok;
        double const d = input.toDouble(&ok);

        if (ok && d >= bottom() && d <= top())
            return QValidator::Acceptable;
        else
            return QValidator::Invalid;
    }
};
