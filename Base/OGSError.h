/**
 * \file OGSError.h
 * KR Initial implementation
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
	static void box(QString e);
	static void box(QString e, QString t);

protected:
	OGSError();
	~OGSError();

};

#endif //OGSERROR_H
