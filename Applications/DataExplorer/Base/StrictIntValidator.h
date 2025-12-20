// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <QIntValidator>

/**
 * \brief A validator for an input field which only accepts integers.
 * Source code adapted from Qt developer faq:
 * https://web.archive.org/web/20100702180529/http://developer.qt.nokia.com/faq/answer/i_can_still_insert_numbers_outside_the_range_specified_with_a_qdoublevalida
 */
class StrictIntValidator : public QIntValidator
{
public:
    StrictIntValidator(int min, int max, QObject* parent = nullptr)
        : QIntValidator(min, max, parent)
    {
    }

    QValidator::State validate(QString& input, int& pos) const override
    {
        if (input.isEmpty())
            return Intermediate;

        if (QIntValidator::validate(input, pos) != Acceptable)
            return Invalid;
        return Acceptable;
    }
};
