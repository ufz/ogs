/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-07
 * \brief  Definition of the MemWatch class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MEMWATCH_H_
#define MEMWATCH_H_

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
    unsigned long _vmem_size = 0;
    unsigned long _rmem_size = 0;
    unsigned long _smem_size = 0;
    unsigned long _cmem_size = 0;
};

}

#endif /* MEMWATCH_H_ */
