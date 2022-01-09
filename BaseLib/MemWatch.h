/**
 * \file
 * \author Thomas Fischer
 * \date   2012-05-07
 * \brief  Definition of the MemWatch class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
