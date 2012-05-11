# OpenGeoSys 6 #

## Intro Base ##

### Logging with logog ###

For details how to use logog see the [OGS devguide](http://ufz.github.com/devguide/logging/) and the [logog documentation](http://johnwbyrd.github.com/logog/).

[logog](http://johnwbyrd.github.com/logog/) is integrated as a [git-subtree](https://github.com/apenwarr/git-subtree) and can be updated with (executed in the sources root):

	git-subtree pull -P BaseLib/logog --squash https://github.com/johnwbyrd/logog.git

It was initially integrated with:

	git-subtree add -P BaseLib/logog --squash https://github.com/johnwbyrd/logog.git master
