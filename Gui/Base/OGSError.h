/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the OGSError class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
