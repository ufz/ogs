// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DemangleTypeInfo.h"

#include <memory>

#ifdef _MSC_VER
// clang-format off
#include <windows.h>
#include <dbghelp.h>
// clang-format on
#else
#include <cxxabi.h>
#endif

namespace BaseLib
{
std::string demangle(const char* mangled_name)
{
#ifdef _MSC_VER
    // MSVC demangling using UnDecorateSymbolName
    char demangled_buffer[1024] = {0};
    if (UnDecorateSymbolName(mangled_name, demangled_buffer,
                             sizeof(demangled_buffer), UNDNAME_COMPLETE))
    {
        return demangled_buffer;
    }
#else
    // GCC/Clang demangling using __cxa_demangle
    int status = 0;
    std::unique_ptr<char, decltype(&std::free)> demangled(
        abi::__cxa_demangle(mangled_name, nullptr, nullptr, &status),
        &std::free);

    if (status == 0 && demangled)
    {
        return demangled.get();
    }
#endif

    // Fallback to mangled name if demangling fails
    return mangled_name;
}

}  // namespace BaseLib
