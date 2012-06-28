/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file DateTools.cpp
 *
 * Created on 2010-06-16 by Karsten Rink
 */

#include "DateTools.h"
#include <cmath>
#include <cstdlib>

namespace BaseLib {

double date2double(int y, int m, int d)
{
	if ( (y<1000 || y>9999) || (m<1 || m>12) || (d<1 || d>31) )
	{
		std::cout << "Error: date2double() -- input not in expected format." << std::endl;
		return 0;
	}

	int ddate=0;
	if (y<1900) y+=1900;
	ddate = y*10000;
	ddate += (m*100);
	ddate += d;

	return ddate;
}

std::string date2string(double ddate)
{
	if (ddate<10000101 || ddate>99991231)
	{
		std::cout << "Error: date2String() -- input not in expected format." << std::endl;
		return "0.0.0000";
	}

	int rest (static_cast<int>(ddate));
	int y = static_cast<int>(floor(rest/10000.0));
	rest = rest % (y*10000);
	int m = static_cast<int>(floor(rest/100.0));
	if (m<1 || m>12) std::cout << "Warning: date2String() -- month not in [1:12]" << std::endl;
	rest = rest % (m*100);
	int d = rest;
	if (d<1 || d>31) std::cout << "Warning: date2String() -- day not in [1:31]" << std::endl;

	std::string day = number2str(d);
	if (d<10) day = "0" + day;
	std::string month = number2str(m);
	if (m<10) month = "0" + month;
	std::string s =  number2str(y) + "-" + month + "-" + day;
	return s;
}

double strDate2double(const std::string &s)
{
	size_t sep ( s.find(".",0) );
	int d ( atoi(s.substr(0, sep).c_str()) );
	size_t sep2 ( s.find(".", sep+1) );
	int m ( atoi(s.substr(sep+1,sep2-(sep+1)).c_str()) );
	int y ( atoi(s.substr(sep2+1, s.length()-(sep2+1)).c_str()) );
	return date2double(y, m, d);
}

double xmlDate2double(const std::string &s)
{
	if (s.length() == 10)
	{
		int d = atoi(s.substr(8,2).c_str());
		if (d<1 || d>31) std::cout << "Warning: xmlDate2double() -- day not in [1:31]" << std::endl;
		int m = atoi(s.substr(5,2).c_str());
		if (m<1 || m>12) std::cout << "Warning: xmlDate2double() -- month not in [1:12]" << std::endl;
		int y = atoi(s.substr(0,4).c_str());
		return date2double(y, m, d);
	}
	return 0;
}

} // end namespace BaseLib
