// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <typeinfo>

namespace BaseLib
{
/// Non-template helper function for demangling mangled type names.
/// Converts a C++ ABI mangled name into a human-readable string.
/// Uses C++ ABI demangling (GCC/Clang) or UnDecorateSymbolName (MSVC).
/// Falls back to the mangled name if demangling fails.
/// See https://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html for
/// GCC/Clang details and
/// https://learn.microsoft.com/en-us/windows/win32/api/dbghelp/nf-dbghelp-undecoratesymbolname
/// for MSVC.
std::string demangle(const char* mangled_name);

/// Get a human-readable type name for error messages.
/// Uses C++ ABI demangling (GCC/Clang) or UnDecorateSymbolName (MSVC) to
/// convert RTTI symbols into readable type names. Falls back to typeid().name()
/// if demangling fails.
///
/// Template overload for explicit type specification:
///   typeToString<SomeType>()
template <typename T>
std::string typeToString()
{
    return demangle(typeid(T).name());
}

/// Template overload for type deduction from values:
///   typeToString(some_value)
/// For non-polymorphic types this equals typeToString<decltype(some_value)>()
/// up to cv-qualifiers and references. For polymorphic types the dynamic type
/// is returned instead, e.g.
/// typeToString(*base_ptr) yields the derived type name; base_ptr must not be
/// null.
template <typename T>
std::string typeToString(T&& value)
{
    return demangle(typeid(value).name());
}
}  // namespace BaseLib
