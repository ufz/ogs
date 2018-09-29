# Third-party libraries #

## logog ##

For details how to use logog see the [OGS devguide](http://devguide.opengeosys.com/logging/) and the [logog documentation](http://johnwbyrd.github.com/logog/).

[logog](http://johnwbyrd.github.com/logog/) is integrated as a [git-subtree](https://github.com/apenwarr/git-subtree) and can be updated with (executed in the sources root):

    git-subtree pull -P ThirdParty/logog --squash https://github.com/johnwbyrd/logog.git

It was initially integrated with:

    git-subtree add -P ThirdParty/logog --squash https://github.com/johnwbyrd/logog.git master

## RapidXML ##

Is used for XML-IO. Is integrated directly.

## tclap ##

Command line option parser. Is integrated directly.

## zlib ##

Compression algorithms. Is integrated directly.

## gtest ##

Google testing framework for unit tests. Is integrated directly.

## autocheck ##

Is integrated as a submodule.

## nlohmann/json ##

Header only c++ json parser. Directly integrated because the benchmarks and
tests are too large.
