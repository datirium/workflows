#################################################################
# Dockerfile
#
# Software:         STAR
# Software Version: 2.5.1b
# Description:      STAR image for BioWardrobe
# Website:          https://github.com/alexdobin/STAR
# Provides:         STAR
# Base Image:       biowardrobe2/scidap:v0.0.3
# Build Cmd:        docker build --rm -t biowardrobe2/star:v2.5.1b -f star-2.5.1b-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/star:v2.5.1b
# Run Cmd:          docker run --rm -ti biowardrobe2/star:v2.5.1b STAR
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.3
MAINTAINER Andrey V Kartashov "porter@porter.st"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Install STAR

ENV VERSION 2.5.1b
ENV NAME STAR
ENV URL "https://github.com/alexdobin/${NAME}/archive/${VERSION}.tar.gz"

RUN wget -q -O - $URL | tar -zxv && \
    cd ${NAME}-${VERSION}/source && \
    make -j 4 && \
    cd .. && \
    cp ./bin/Linux_x86_64_static/${NAME} /usr/local/bin/ && \
    cd .. && \
    strip /usr/local/bin/${NAME}; true && \
    rm -rf ./${NAME}-${VERSION}/
