/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#if !defined(_WIN32) && !defined(__APPLE__)
#include <cfenv>

class DisableFPE
{
public:
    DisableFPE()
    {
        fegetenv(&fe_env);     // Store floating-point exception handling.
        fesetenv(FE_DFL_ENV);  // Set default environment effectively disabling
                               // exceptions.
    }

    ~DisableFPE()
    {
        fesetenv(&fe_env);  // Restore floating-point exception handling.
    }

private:
    fenv_t fe_env;
};

#else

// No-op for Windows and Apple
class DisableFPE
{
};

#endif
