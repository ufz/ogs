/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file DateTools.h
 *
 * Created on 2010-01-22 by Karsten Rink
 */

#ifndef DATETOOLS_H
#define DATETOOLS_H

#include "StringTools.h"
#include <string>

namespace BaseLib {

/**
 * Converts three integers representing a date into a double.
 * Note: It is not really checked if the date actually makes sense.
 */
int date2int(int y, int m, int d);

/**
 * Converts an integer to a string date "dd.mm.yyyy"
 * Note: (Almost) no checks are performed if the int makes sense as a date.
 */
std::string int2date(int date);

/**
 * Converts a double representing a date into a string.
 * Note: It is not really checked if the date actually makes sense.
 * \param ddate Number containing date in double format yyyymmdd
 * \return A string containing the date in format "dd.mm.yyyy".
 */
std::string date2string(double ddate);

/**
 * Converts a string containing a date into a double.
 * Note: It is not really checked if the date actually makes sense.
 * \param s String containing the date, the expected format is "dd.mm.yyyy".
 * \return A number representing the date as dd.mm.yyyy.
 */
int strDate2int(const std::string &s);

/**
 * Converts a string containing a date into a double.
 * Note: It is not really checked if the date actually makes sense.
 * \param s String containing the date, the expected format is conform to the xml date type, i.e. "yyyy-mm-dd".
 * \return A number representing the date as yyyymmdd.
 */
int xmlDate2int(const std::string &s);

} // namespace BaseLib

#endif //DATETOOLS_H
