Defines for a *dynamic* secondary variable which of its input data are provided
by which function.

Said function must be defined in OGS as a NumLib::NamedFunction and must be available to
the current process.
In order to find out which named functions are available, you can search where
the method NumLib::NamedFunctionProvider::getNamedFunctions() is implemented.
