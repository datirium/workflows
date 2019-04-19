#################################################################
# Dockerfile
#
# Software:         geep
# Software Version: 0.0.4
# Description:      RPKM calculation tool (includes bamtools)
# Website:          https://github.com/pezmaster31/bamtools, https://github.com/Barski-lab/geep/tree/master
# Provides:         geep
# Base Image:       biowardrobe2/scidap:v0.0.3
# Build Cmd:        docker build --rm -t biowardrobe2/geep:v0.0.4 -f geep-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/geep:v0.0.4
# Run Cmd:          docker run --rm -ti biowardrobe2/geep:v0.0.4
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.3
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing geep

ENV VERSION 0.0.4
ENV NAME geep
ENV URL "https://github.com/Barski-lab/${NAME}/archive/v${VERSION}.tar.gz"
ENV LD_LIBRARY_PATH "/usr/local/lib/bamtools:${LD_LIBRARY_PATH}"


RUN apt-get update && \
    apt-get install libboost-all-dev -y && \
    wget -q -O - $URL | tar -zxv && \
    cd ${NAME}-${VERSION} && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make install && \
    cd ../.. && \
    rm -rf ./${NAME}-${VERSION}/ && \
    strip /usr/local/bin/*; true
