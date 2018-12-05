FROM ubuntu:17.10

RUN apt-get update && apt-get install -y curl \
  && curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash \
  && apt-get update && apt-get -y install \
    ccache \
    g++-6 \
    git git-lfs \
    python3-pip \
    sudo \
    unzip

ENV CC=gcc-6
ENV CXX=g++-6

RUN python3 -m pip install --upgrade pip \
  && python3 -m pip install cmake conan>=1.10.0

# Ninja
RUN curl -L -o ninja-linux.zip https://github.com/ninja-build/ninja/releases/download/v1.8.2/ninja-linux.zip \
  && unzip ninja-linux.zip \
  && mv ninja /usr/local/bin/ninja \
  && rm ninja-linux.zip

 # Add user jenkins to the image
RUN adduser --uid 500 --disabled-password --gecos "" jenkins \
  # Add user jenkins to sudoers with NOPASSWD
  && echo "jenkins ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
  # Set password for the jenkins user (you may want to alter this).
  && echo "jenkins:jenkins" | chpasswd

USER jenkins
ENV CCACHE_DIR=/home/jenkins/cache/ccache
RUN mkdir -p $CCACHE_DIR
ENV CCACHE_MAXSIZE=15G
ENV CCACHE_SLOPPINESS=pch_defines,time_macros
WORKDIR /home/jenkins
RUN conan user
