// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
