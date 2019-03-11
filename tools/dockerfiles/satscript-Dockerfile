#################################################################
# Dockerfile
#
# Software:         SatScript
# Software Version: v0.0.2
# Description:      SatScript image for BioWardrobe plugin to plot Saturation plots
# Website:          https://github.com/michael-kotliar/SatScript
# Provides:         SatScript
#                   bedops
#                   bedmap
#                   bc
# Base Image:       biowardrobe2/macs2:v2.1.1
# Build Cmd:        docker build --rm -t biowardrobe2/satscript:v0.0.2 -f satscript-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/satscript:v0.0.2
# Run Cmd:          docker run --rm -ti biowardrobe2/satscript:v0.0.2 SatScript --help
#################################################################
# v0.0.2
# - update SatScript v0.0.2 allows to calculate FRIP
# - add Bedops
#################################################################


### Base Image
FROM biowardrobe2/macs2:v2.1.1
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive


################## BEGIN INSTALLATION ######################

WORKDIR /tmp

ENV MPLBACKEND "AGG"

ENV SATSCRIPT_VERSION 0.0.2
ENV SATSCRIPT_URL "https://github.com/michael-kotliar/SatScript/archive/v${SATSCRIPT_VERSION}.tar.gz"

ENV BEDOPS_VERSION 2.4.35
ENV BEDOPS_URL "https://github.com/bedops/bedops/archive/v${BEDOPS_VERSION}.tar.gz"


### Installing required packages

RUN apt-get clean all &&\
    apt-get update &&\
    apt-get install -y libc6-dev build-essential tcsh devscripts debhelper bc &&\

### Installing bedops
    wget -q -O - $BEDOPS_URL | tar -zxv --strip-components=1 &&\
    make -j `nproc` &&\
    make install BINDIR=/usr/local/bin/ &&\

### Installing SatScript
    wget -q -O - $SATSCRIPT_URL | tar -zxv --strip-components=1 &&\
    pip install . &&\

### Cleaning
    apt-get clean &&\
    apt-get purge &&\
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* &&\
    strip /usr/local/bin/*; true