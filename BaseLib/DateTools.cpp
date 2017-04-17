/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-16
 * \brief  Implementation of date helper functions.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DateTools.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <sstream>

#include <logog/include/logog.hpp>

namespace BaseLib
{
int date2int(int y, int m, int d)
{
    if ( (y < 1000 || y > 9999) || (m < 1 || m > 12) || (d < 1 || d > 31) )
    {
        WARN("date2int(): Input not in expected range.");
        return 0;
    }

    int ddate(0);
    ddate = y * 10000;
    ddate += (m * 100);
    ddate += d;

    return ddate;
}

std::string int2date(int date)
{
    if (date > 10000000 && date < 22000000)
    {
        auto y = static_cast<int>(std::floor(date / 10000.0));
        auto m = static_cast<int>(std::floor((date - (y * 10000)) / 100.0));
        int d = date - (y * 10000) - (m * 100);
        std::stringstream ss;
        if (d < 10)
            ss << "0";
        ss << d << ".";
        if (m < 10)
            ss << "0";
        ss << m << "." << y;
        return ss.str();
    }
    return "";
}

std::string date2string(double ddate)
{
    if (ddate < 10000101 || ddate > 99991231)
    {
        WARN("date2String(): Input not in expected format.");
        return "0.0.0000";
    }

    auto rest(static_cast<int>(ddate));
    auto y = static_cast<int>(std::floor(rest / 10000.0));
    rest = rest % (y * 10000);
    auto m = static_cast<int>(std::floor(rest / 100.0));
    if (m < 1 || m > 12)
        WARN("date2String(): month not in [1:12].");
    rest = rest % (m * 100);
    int d = rest;
    if (d < 1 || d > 31)
        WARN("date2String(): day not in [1:31].");

    std::string day = std::to_string(d);
    if (d < 10)
        day = "0" + day;
    std::string month = std::to_string(m);
    if (m < 10)
        month = "0" + month;
    std::string s =  std::to_string(y) + "-" + month + "-" + day;
    return s;
}

int strDate2int(const std::string &s)
{
    std::string str(s);
    if (s.length() > 10)
        str = s.substr(0,10);
    std::size_t sep ( str.find('.',0) );
    int d ( atoi(str.substr(0, sep).c_str()) );
    std::size_t sep2 ( str.find('.', sep + 1) );
    int m ( atoi(str.substr(sep + 1,sep2 - (sep + 1)).c_str()) );
    int y ( atoi(str.substr(sep2 + 1, s.length() - (sep2 + 1)).c_str()) );
    return date2int(y, m, d);
}

int xmlDate2int(const std::string &s)
{
    if (s.length() == 10)
    {
        int d = atoi(s.substr(8,2).c_str());
        if (d < 1 || d > 31)
            WARN("xmlDate2double(): day not in [1:31].");
        int m = atoi(s.substr(5,2).c_str());
        if (m < 1 || m > 12)
            WARN("xmlDate2double(): month not in [1:12].");
        int y = atoi(s.substr(0,4).c_str());
        return date2int(y, m, d);
    }
    return 0;
}

std::string formatDate(
    std::chrono::time_point<std::chrono::system_clock> const& time)
{
    auto const time_t = std::chrono::system_clock::to_time_t(time);
    char time_str[100];
    if (std::strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S%z",
                      std::localtime(&time_t))) {
        return time_str;
    }
    return "FAILED FORMATTING THE GIVEN TIME POINT.";
}
} // end namespace BaseLib
