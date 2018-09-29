Determines if stdout should be flushed each time before and after passing control to the Python script.
This might slow down operation, but ensures the right order of messages printed
to stdout from Python and C++.
