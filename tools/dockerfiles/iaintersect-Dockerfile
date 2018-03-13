#################################################################
# Dockerfile
#
# Software:         iaintersect
# Software Version: 0.0.2
# Description:      iaintersect image for BioWardrobe
# Website:          https://github.com/Barski-lab/iaintersect
# Provides:         iaintersect
# Base Image:       biowardrobe2/iaintersect:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/iaintersect:v0.0.2 -f iaintersect-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/iaintersect:v0.0.2
# Run Cmd:          docker run --rm -ti biowardrobe2/iaintersect:v0.0.2 iaintersect
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing iaintersect

ENV VERSION 0.0.2
ENV NAME iaintersect
ENV URL "https://github.com/Barski-lab/${NAME}/archive/v${VERSION}.tar.gz"

RUN apt-get update && \
    apt-get install qt5-default -y && \
    apt-get install libboost-dev -y && \
    wget -q -O - $URL | tar -zxv && \
    cd ${NAME}-${VERSION} && \
    qmake && \
    make && \
    cp ./${NAME} /usr/local/bin/ && \
    strip /usr/local/bin/*; true && \
    cd .. && \
    rm -rf ./${NAME}-${VERSION}/