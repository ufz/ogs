// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <vector>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
/// Manages a pool of independent IPhreeqc instances for thread-safe parallel
/// chemistry calculations.
///
/// Each OpenMP thread uses its own PHREEQC instance to avoid thread-safety
/// issues. All instances are initialized with the same thermodynamic database.
class PhreeqcInstancePool
{
public:
    /// Creates a pool of PHREEQC instances.
    /// @param database Path to the thermodynamic database file.
    /// @param num_instances Number of instances to create (one per thread).
    PhreeqcInstancePool(std::string const& database, int num_instances);

    ~PhreeqcInstancePool();

    // Disable copy
    PhreeqcInstancePool(PhreeqcInstancePool const&) = delete;
    PhreeqcInstancePool& operator=(PhreeqcInstancePool const&) = delete;

    /// Returns the PHREEQC instance ID for the given thread.
    /// @param thread_id The OpenMP thread ID (0-based).
    /// @return The IPhreeqc instance ID to use for this thread.
    int getInstanceForThread(int thread_id) const;

    /// Returns the number of instances in the pool.
    int size() const { return static_cast<int>(instance_ids_.size()); }

    /// Creates a new IPhreeqc instance and loads the thermodynamic database.
    /// Returns the instance id. Fatals on failure.
    static int createInstance(std::string const& database);

    /// Configures an IPhreeqc instance for string-based I/O: disables file
    /// output and enables in-memory buffers for selected output and dump.
    static void configureForStringIO(int id);

private:
    std::vector<int> instance_ids_;
};

}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
