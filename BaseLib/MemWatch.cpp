/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-07
 * \brief  Implementation of the MemWatch class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MemWatch.h"

#if !defined(_WIN32) && !defined(__APPLE__) && !defined(__MINGW32__)
#include <stdio.h>
#include <unistd.h>

#include <fstream>
#include <sstream>
#include <string>
#endif

namespace BaseLib
{
MemWatch::MemWatch()
{
    updateMemUsage();
}

unsigned MemWatch::updateMemUsage()
{
#if !defined(_WIN32) && !defined(__APPLE__) && !defined(__MINGW32__)
    std::string fname("/proc/");
    std::stringstream str_pid;
    str_pid << static_cast<unsigned>(getpid());
    fname += str_pid.str();
    fname += "/statm";
    unsigned pages;

    std::ifstream in(fname.c_str(), std::ios::in);
    if (!in.is_open())
    {
        perror("open");
        return 1;
    }

    in >> pages;
    vmem_size_ = static_cast<unsigned long>(pages) *
                 static_cast<unsigned long>(getpagesize());
    in.close();
#endif

    return 0;
}

unsigned long MemWatch::getVirtMemUsage()
{
    updateMemUsage();
    return vmem_size_;
}

}  // end namespace BaseLib
