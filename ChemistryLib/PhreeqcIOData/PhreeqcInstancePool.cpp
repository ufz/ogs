// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PhreeqcInstancePool.h"

#include <IPhreeqc.h>

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
int PhreeqcInstancePool::createInstance(std::string const& database)
{
    int const id = CreateIPhreeqc();
    if (id < 0)
    {
        OGS_FATAL(
            "Failed to create PHREEQC instance (IPQ_RESULT error code = {}).",
            id);
    }
    if (LoadDatabase(id, database.c_str()) != IPQ_OK)
    {
        OutputErrorString(id);
        OGS_FATAL(
            "Failed to load thermodynamic database '{:s}' for PHREEQC "
            "instance {}.",
            database, id);
    }
    return id;
}

void PhreeqcInstancePool::configureForStringIO(int const id)
{
    SetSelectedOutputFileOn(id, 0);
    if (SetSelectedOutputStringOn(id, 1) != IPQ_OK)
    {
        OGS_FATAL("Failed to enable string output for PHREEQC instance {}.",
                  id);
    }
    SetDumpFileOn(id, 0);
    if (SetDumpStringOn(id, 1) != IPQ_OK)
    {
        OGS_FATAL("Failed to enable dump string for PHREEQC instance {}.", id);
    }
}

PhreeqcInstancePool::PhreeqcInstancePool(std::string const& database,
                                         int const num_instances)
{
    if (num_instances < 1)
    {
        OGS_FATAL("PhreeqcInstancePool requires at least 1 instance, got {}.",
                  num_instances);
    }

    instance_ids_.reserve(num_instances);
    for (int i = 0; i < num_instances; ++i)
    {
        int const id = createInstance(database);
        configureForStringIO(id);
        instance_ids_.push_back(id);
    }

    INFO("Created {} PHREEQC instances for parallel chemistry execution.",
         num_instances);
}

PhreeqcInstancePool::~PhreeqcInstancePool()
{
    for (int const id : instance_ids_)
    {
        DestroyIPhreeqc(id);
    }
}

int PhreeqcInstancePool::getInstanceForThread(int thread_id) const
{
    if (thread_id < 0 || thread_id >= static_cast<int>(instance_ids_.size()))
    {
        OGS_FATAL(
            "Thread ID {} is out of range [0, {}). This indicates a bug in "
            "the parallel chemistry implementation.",
            thread_id, instance_ids_.size());
    }
    return instance_ids_[thread_id];
}

}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
