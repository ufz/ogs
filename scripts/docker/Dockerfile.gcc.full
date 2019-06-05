FROM ubuntu:18.04

# Generated with ogs-container-maker 1.2.0

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        curl \
        tar \
        wget && \
    rm -rf /var/lib/apt/lists/*

# GNU compiler
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        g++ \
        gcc && \
    rm -rf /var/lib/apt/lists/*

# OGS base building block
# Python
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python \
        python-dev \
        python3 \
        python3-dev && \
    rm -rf /var/lib/apt/lists/*
# pip
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3-pip \
        python3-setuptools \
        python3-wheel && \
    rm -rf /var/lib/apt/lists/*
RUN pip3 install virtualenv
# CMake version 3.13.4
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        wget && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://cmake.org/files/v3.13/cmake-3.13.4-Linux-x86_64.sh && \
    /bin/sh /var/tmp/cmake-3.13.4-Linux-x86_64.sh --prefix=/usr/local --skip-license && \
    rm -rf /var/tmp/cmake-3.13.4-Linux-x86_64.sh
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends software-properties-common && \
    apt-add-repository ppa:git-core/ppa -y && \
    apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        git \
        git-lfs \
        make \
        ninja-build && \
    rm -rf /var/lib/apt/lists/*
RUN apt-get update && \
    apt-get install -y dirmngr --install-recommends && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 6B05F25D762E3157 && \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && \
    git lfs install && \
    mkdir -p /apps /scratch /lustre /work /projects

# Package manager Conan building block
# pip
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3-pip \
        python3-setuptools \
        python3-wheel && \
    rm -rf /var/lib/apt/lists/*
RUN pip3 install conan==1.16.0
RUN mkdir -p /opt/conan && \
    chmod 777 /opt/conan
ENV CONAN_USER_HOME=/opt/conan
LABEL org.opengeosys.pm=conan \
    org.opengeosys.pm.conan.version=1.16.0
LABEL org.opengeosys.pm.conan.user_home=/opt/conan

# cppcheck version 1.87

RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://github.com/danmar/cppcheck/archive/1.87.tar.gz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/1.87.tar.gz -C /var/tmp -z && \
    mkdir -p /var/tmp/build && cd /var/tmp/build && cmake -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_INSTALL_PREFIX=/usr/local/cppcheck -DCMAKE_BUILD_TYPE=RELEASE /var/tmp/cppcheck-1.87 && \
    cmake --build /var/tmp/build --target install -- -j$(nproc) && \
    rm -rf /var/tmp/1.87.tar.gz /var/tmp/build /var/tmp/cppcheck-1.87
ENV PATH=/usr/local/cppcheck/bin:$PATH

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        doxygen \
        graphviz \
        texlive-base && \
    rm -rf /var/lib/apt/lists/*

# pip
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3-pip \
        python3-setuptools \
        python3-wheel && \
    rm -rf /var/lib/apt/lists/*
RUN pip3 install gcovr

# Package manager Conan building block
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        ccache && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir -p /opt/cache && chmod 777 /opt/cache
ENV CCACHE_DIR=/opt/cache \
    CCACHE_MAXSIZE=15G \
    CCACHE_SLOPPINESS=pch_defines,time_macros
LABEL ccache.dir=/opt/cache \
    ccache.size=15G

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        sudo && \
    rm -rf /var/lib/apt/lists/*

# Jenkins node
RUN groupadd --gid 1001 jenkins && \
    adduser --uid 500 --gid 1001 --disabled-password --gecos "" jenkins && \
    echo "jenkins ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers && \
    echo "jenkins:jenkins" | chpasswd
USER jenkins
WORKDIR /home/jenkins
