/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
