/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-07
 * \brief  Implementation of the MemWatch class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MemWatch.h"

#if !defined(WIN32) && !defined(__APPLE__) && !defined(__MINGW32__)
#include <fstream>
#include <string>
#include <sstream>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace BaseLib {

MemWatch::MemWatch ()
{
        updateMemUsage ();
}

unsigned MemWatch::updateMemUsage ()
{
#if !defined(WIN32) && !defined(__APPLE__) && !defined(__MINGW32__)
        std::string fname ("/proc/");
        std::stringstream str_pid;
        str_pid << static_cast<unsigned> (getpid());
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
        _vmem_size = static_cast<unsigned long>(pages) *
            static_cast<unsigned long>(getpagesize());
        in >> pages;
        _rmem_size =static_cast<unsigned long>(pages) *
            static_cast<unsigned long>(getpagesize());
        in >> pages;
        _smem_size = static_cast<unsigned long>(pages) *
            static_cast<unsigned long>(getpagesize());
        in >> pages;
        _cmem_size = static_cast<unsigned long>(pages) *
            static_cast<unsigned long>(getpagesize());
        in.close ();
#endif

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

