// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>
#include <string>
#include <string_view>

namespace ProcessLib
{
/// Interprets the \c coupling_scheme value read from a process configuration.
///
/// \param coupling_scheme the configured value, or nothing if the process
/// configuration does not contain the tag.
/// \param process_name the process type, used in the error message only.
///
/// \return true for the monolithic scheme, which is also the scheme used if
/// \c coupling_scheme is unset. The deprecated "monolithic_scheme" is accepted
/// as "monolithic" with a warning.
///
/// Any value other than "monolithic" or "staggered" is fatal (\c OGS_FATAL):
/// a misspelled scheme must not silently select the monolithic one.
bool useMonolithicScheme(std::optional<std::string> const& coupling_scheme,
                         std::string_view const process_name);
}  // namespace ProcessLib
