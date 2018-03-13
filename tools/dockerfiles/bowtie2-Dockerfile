#################################################################
# Dockerfile
#
# Software:         bowtie2
# Software Version: 2.3.0
# Description:      Bowtie2 image for BioWardrobe
# Website:          http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# Provides:         bowtie2
# Base Image:       biowardrobe2/scidap:v0.0.2
# Build Cmd:        docker build --rm -t biowardrobe2/bowtie2:v2.3.0 -f bowtie2-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/bowtie2:v2.3.0
# Run Cmd:          docker run --rm -ti biowardrobe2/bowtie2:v2.3.0 bowtie2
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.2
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Installing bowtie2

ENV VERSION 2.3.0
ENV NAME bowtie2
ENV URL "https://github.com/BenLangmead/bowtie2/archive/v${VERSION}.tar.gz"

### First installing libtbb, then Bowtie2

RUN apt-get update && apt-get install -y libtbb-dev && \
    wget -q -O - $URL | tar -zxv && \
    cd ${NAME}-${VERSION} && \
    make -j 4 && \
    cd .. && \
    cp ./${NAME}-${VERSION}/${NAME} /usr/local/bin/ && \
    cp ./${NAME}-${VERSION}/${NAME}-* /usr/local/bin/ && \
    strip /usr/local/bin/*; true && \
    rm -rf ./${NAME}-${VERSION}/
