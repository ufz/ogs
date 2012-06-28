/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file FileTools.h
 *
 * Created on 2010-04-26 by Lars Bilke
 *
 */


#ifndef FILETOOLS_H
#define FILETOOLS_H

// ** INCLUDES **
#include <sys/stat.h>

namespace BaseLib {
/**
 * Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 */
static bool IsFileExisting(std::string strFilename)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;

	// Attempt to get the file attributes
	intStat = stat(strFilename.c_str(),&stFileInfo);

	if(intStat == 0)
	{
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	}
	else
	{
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}

	return(blnReturn);
}

} // end namespace BaseLib

#endif // FILETOOLS_H
