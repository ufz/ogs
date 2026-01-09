// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

// A macro preventing inlining of a function, to guide compiler optimizations
#ifdef _MSC_VER
#define OGS_NO_INLINE __declspec(noinline)
#else
#define OGS_NO_INLINE __attribute__((noinline))
#endif
