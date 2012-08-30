/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file OGSError.h
 *
 * Created on by Karsten Rink
 */

#ifndef OGSERROR_H
#define OGSERROR_H

class QString;

/**
 * \brief Manages error messages via message boxes
 */
class OGSError
{
public:
	static void box(const QString &e);
	static void box(const QString &e, const QString &t);

protected:
	OGSError();
	~OGSError();
};

#endif //OGSERROR_H
