// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#if defined(__GNUC__) && __GNUC__ >= 13
#define OGS_NO_DANGLING [[gnu::no_dangling]]
#else
#define OGS_NO_DANGLING
#endif
