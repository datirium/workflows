#################################################################
# Dockerfile
#
# Software:         atdp
# Software Version: 0.0.1
# Description:      atdp image for BioWardrobe
# Website:          https://github.com/Barski-lab/atdp
# Provides:         atdp, refgene-sort
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/atdp:v0.0.1 -f atdp-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/atdp:v0.0.1
# Run Cmd:          docker run --rm -ti biowardrobe2/atdp:v0.0.1 atdp
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp



ENV VERSION_ATDP 0.0.1
ENV NAME_ATDP atdp
ENV URL_ATDP "https://github.com/Barski-lab/${NAME_ATDP}/archive/v${VERSION_ATDP}.tar.gz"

ENV VERSION_SORT 0.0.1
ENV NAME_SORT refgene-sort
ENV URL_SORT "https://github.com/Barski-lab/${NAME_SORT}/archive/v${VERSION_SORT}.tar.gz"

### Installing atdp

RUN apt-get update && \
    apt-get install qt5-default -y && \
    apt-get install libboost-dev -y && \
    wget -q -O - $URL_ATDP | tar -zxv && \
    cd ${NAME_ATDP}-${VERSION_ATDP} && \
    qmake && \
    make && \
    cp ./${NAME_ATDP} /usr/local/bin/ && \
    strip /usr/local/bin/*; true && \
    cd .. && \
    rm -rf ./${NAME_ATDP}-${VERSION_ATDP}/ && \

### Installing refgene-sort

    wget -q -O - $URL_SORT | tar -zxv && \
    cd ${NAME_SORT}-${VERSION_SORT} && \
    pip install setuptools && \
    python setup.py install && \
    cd .. && \
    rm -rf ./${NAME_SORT}-${VERSION_SORT}
