#################################################################
# Dockerfile
#
# Software:         fastx_toolkit
# Software Version: 0.0.14
# Description:      Fastx toolkit for BioWardrobe
# Website:          http://hannonlab.cshl.edu/fastx_toolkit/
# Provides:         fastx_toolkit
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/fastx_toolkit:v0.0.14 -f fastx-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/fastx_toolkit:v0.0.14
# Run Cmd:          docker run --rm -ti biowardrobe2/fastx_toolkit:v0.0.14 fastx_quality_stats -h
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Andrey V Kartashov "porter@porter.st"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Install fastx_toolkit

ENV VERSION 0.0.14
ENV NAME fastx_toolkit
ENV URL "https://github.com/agordon/fastx_toolkit/releases/download/${VERSION}/${NAME}-${VERSION}.tar.bz2"

RUN apt-get clean all &&\
    apt-get update &&\
    apt-get install -y  \
        libgtextutils-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget -q -O - $URL | tar -jxv && \
    cd ${NAME}-${VERSION} && \
    ./configure --prefix=/usr/local/ && \
    make -j 4 && \
    make install && \
    cd .. && \
    rm -rf ./${NAME}-${VERSION}/ && \
    strip /usr/local/bin/*; true
