// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <QDateTime>
#include <QString>
#include <string>

/// Converts string into QDateTime-format
inline QDateTime getDateTime(QString const& stringDate)
{
    if (stringDate.length() == 10)
        return QDateTime::fromString(stringDate, "dd.MM.yyyy");

    if (stringDate.length() == 19)
        return QDateTime::fromString(stringDate, "dd.MM.yyyy.HH.mm.ss");

    return QDateTime();
}

/// Converts string into QDateTime-format
inline QDateTime getDateTime(std::string const& stringDate)
{
    return getDateTime(QString::fromStdString(stringDate));
}
