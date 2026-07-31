// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#if !defined(_WIN32) && !defined(__APPLE__)
// Some libstdc++ toolchains do not detect the platform's <fenv.h> during
// their build. Force the C++ wrapper to include it in that case.
#if defined(__linux__) && defined(__GNUC__) && !defined(_GLIBCXX_HAVE_FENV_H)
#define _GLIBCXX_HAVE_FENV_H 1
#include <fenv.h>
#undef _GLIBCXX_HAVE_FENV_H
#else
#include <fenv.h>
#endif
#endif
