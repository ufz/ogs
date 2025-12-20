// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
