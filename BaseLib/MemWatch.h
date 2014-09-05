/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-07
 * \brief  Definition of the MemWatch class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MEMWATCH_H_
#define MEMWATCH_H_

#ifndef _MSC_VER

#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

namespace BaseLib {

class MemWatch {
public:
	MemWatch ();
	unsigned long getVirtMemUsage ();
	unsigned long getResMemUsage ();
	unsigned long getShrMemUsage ();
	unsigned long getCodeMemUsage ();

private:
	unsigned updateMemUsage ();
	unsigned long _vmem_size;
	unsigned long _rmem_size;
	unsigned long _smem_size;
	unsigned long _cmem_size;
};

}

#endif // not Windows

#endif /* MEMWATCH_H_ */
