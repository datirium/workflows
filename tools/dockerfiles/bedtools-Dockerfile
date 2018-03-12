#################################################################
# Dockerfile
#
# Software:         bedtools2
# Software Version: 2.26.0
# Description:      bedtools2 image for BioWardrobe
# Website:          http://bedtools.readthedocs.org/
# Provides:         bedtools2
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/bedtools2:v2.26.0 -f bedtools-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/bedtools2:v2.26.0
# Run Cmd:          docker run --rm -ti biowardrobe2/bedtools2:v2.26.0 bedtools
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Andrey V Kartashov "porter@porter.st"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing bedtools2

ENV VERSION 2.26.0
ENV NAME bedtools2
ENV URL "https://github.com/arq5x/bedtools2/releases/download/v${VERSION}/bedtools-${VERSION}.tar.gz"

RUN wget -q -O - $URL | tar -zxv && \
cd ${NAME} && \
make -j 4 && \
cd .. && \
cp ./${NAME}/bin/bedtools /usr/local/bin/ && \
strip /usr/local/bin/*; true && \
rm -rf ./${NAME}/
