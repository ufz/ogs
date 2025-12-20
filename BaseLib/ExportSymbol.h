// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#ifndef OGS_EXPORT_SYMBOL
#if defined(WIN32) || defined(_WIN32)
#define OGS_EXPORT_SYMBOL __declspec(dllexport)
#else
#define OGS_EXPORT_SYMBOL __attribute__((visibility("default")))
#endif
#endif  // defined(OGS_EXPORT_SYMBOL)
