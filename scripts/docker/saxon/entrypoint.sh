#!/bin/bash

java $JAVA_OPTS -jar /usr/src/saxon/Saxon-HE-9.9.1-5.jar -xsl:/usr/src/saxon/ctest-to-junit.xsl "$@"
