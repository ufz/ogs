// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
