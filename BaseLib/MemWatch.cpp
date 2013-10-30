/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-07
 * \brief  Implementation of the MemWatch class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MemWatch.h"

#ifndef _MSC_VER

namespace BaseLib {

MemWatch::MemWatch ()
{
        updateMemUsage ();
}

unsigned MemWatch::updateMemUsage ()
{
        std::string fname ("/proc/");
        std::stringstream str_pid;
        str_pid << (unsigned) getpid();
        fname += str_pid.str();
        fname += "/statm";
        unsigned pages;

        std::ifstream in (fname.c_str(), std::ios::in);
        if (!in.is_open())
        {
            perror( "open" );
            return 1;
        }

        in >> pages;
        _vmem_size = ((unsigned long) pages) * ((unsigned long) getpagesize());
        in >> pages;
        _rmem_size = ((unsigned long) pages) * ((unsigned long) getpagesize());
        in >> pages;
        _smem_size = ((unsigned long) pages) * ((unsigned long) getpagesize());
        in >> pages;
        _cmem_size = ((unsigned long) pages) * ((unsigned long) getpagesize());
        in.close ();
        return 0;
}

unsigned long MemWatch::getVirtMemUsage ()
{
        updateMemUsage ();
        return _vmem_size;
}

unsigned long MemWatch::getResMemUsage () {
        updateMemUsage ();
        return _rmem_size;
}

unsigned long MemWatch::getShrMemUsage () {
        updateMemUsage ();
        return _smem_size;
}

unsigned long MemWatch::getCodeMemUsage () {
        updateMemUsage ();
        return _cmem_size;
}

} // end namespace BaseLib

#endif // _MSC_VER
