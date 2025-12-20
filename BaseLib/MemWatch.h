// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace BaseLib
{
class MemWatch
{
public:
    MemWatch();
    unsigned long getVirtMemUsage();

private:
    unsigned updateMemUsage();
    unsigned long vmem_size_ = 0;
};

}  // namespace BaseLib
