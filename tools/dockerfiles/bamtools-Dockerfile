#################################################################
# Dockerfile
#
# Software:         bamtools
# Software Version: 2.4.1
# Description:      Bamtools image for BioWardrobe
# Website:          https://github.com/pezmaster31/bamtools
# Provides:         bamtools
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/bamtools:v2.4.1 -f bamtools-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/bamtools:v2.4.1
# Run Cmd:          docker run --rm -ti biowardrobe2/bamtools:v2.4.1 bamtools
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing bamtools

ENV VERSION 2.4.1
ENV NAME bamtools
ENV URL "https://github.com/pezmaster31/bamtools/archive/v${VERSION}.tar.gz"

RUN wget -q -O - $URL | tar -zxv && \
    cd ${NAME}-${VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    cd ../.. && \
    cp ./${NAME}-${VERSION}/lib/libbamtools.so.${VERSION} /usr/lib/ && \
    rm -rf ./${NAME}-${VERSION} && \
    strip /usr/local/bin/*; true
