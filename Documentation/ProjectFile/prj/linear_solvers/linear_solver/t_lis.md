Defines an option string that will be passed to solvers of the LIS library.

The default is <tt>-initx_zeros 0</tt>, i.e., not to overwrite the passed initial
guess with zeroes. All other options are appended to this default value.

\attention
LIS can also receive options from the commandline arguments. However, the values
from the prj file are applied directly to the linear solver. If specifying both,
make sure that you check which options take precedence.

Refer to the [documentation of LIS](http://www.ssisc.org/lis/index.en.html) for further details.
