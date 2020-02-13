FROM openjdk:8-jre

WORKDIR /usr/src/saxon
RUN curl -L -O https://repo1.maven.org/maven2/net/sf/saxon/Saxon-HE/9.9.1-5/Saxon-HE-9.9.1-5.jar
COPY entrypoint.sh entrypoint.sh
ENTRYPOINT ["/bin/sh", "/usr/src/saxon/entrypoint.sh"]

RUN curl -L -O https://raw.githubusercontent.com/rpavlik/jenkins-ctest-plugin/master/ctest-to-junit.xsl
