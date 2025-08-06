/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <dlfcn.h>

#include <cassert>
#include <iostream>

void* load(char const* lib)
{
    void* handle = dlopen(lib, RTLD_NOW | RTLD_GLOBAL);

    if (!handle)
    {
        std::cerr << dlerror() << '\n';
        std::exit(EXIT_FAILURE);
    }

    dlerror(); /* Clear any existing error */

    return handle;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "you invoked this command without arguments\n";
        std::exit(EXIT_FAILURE);
    }

    for (int i = 1; i < argc; ++i)
    {
        void* handle = load(argv[i]);
        std::cout << "lib address of " << argv[i] << " is: " << handle
                  << std::endl;
    }
}
