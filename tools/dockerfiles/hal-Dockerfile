#################################################################
# Dockerfile
#
# Software:         HAL
# Software Version: commit 6541a8ebd912f65ed0ae22e96dbb03240b5d069e
# Description:      HAL image for BioWardrobe
# Website:          https://github.com/ComparativeGenomicsToolkit/hal
# Provides:         HAL tools
# Base Image:       biowardrobe2/scidap:v0.0.3
# Build Cmd:        docker build --rm -t biowardrobe2/hal:v0.0.1 -f hal-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/hal:v0.0.1
# Run Cmd:          docker run --rm -ti biowardrobe2/hal:v0.0.1
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.3
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

ENV PATH "$PATH:/usr/local/bin/hdf5/bin:/opt/hal/bin"
ENV PYTHONPATH "/opt:$PYTHONPATH"
ENV h5prefix "-prefix=/usr/local/bin/hdf5"

ENV VERSION_HAL_COMMIT 6541a8ebd912f65ed0ae22e96dbb03240b5d069e
ENV URL_HAL "https://github.com/ComparativeGenomicsToolkit/hal.git"

ENV VERSION_HDF5 1.10
ENV SUB_VERSION_HDF5 2
ENV URL_HDF5 "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${VERSION_HDF5}/hdf5-${VERSION_HDF5}.${SUB_VERSION_HDF5}/src/hdf5-${VERSION_HDF5}.${SUB_VERSION_HDF5}.tar.gz"

ENV VERSION_SONLIB fe6cf4527b3716b983f117a37b31c6e2c8ac6cbe
ENV URL_SONLIB "https://github.com/ComparativeGenomicsToolkit/sonLib.git"

### Install required packages

RUN apt-get clean all &&\
    apt-get update &&\

    cd /tmp &&\
    wget -q -O - ${URL_HDF5} | tar -zxv --strip-components=1 &&\
    ./configure --enable-cxx --prefix /usr/local/bin/hdf5 &&\
    make && make install &&\
    rm -rf /tmp/* &&\

    cd /opt &&\
    git clone ${URL_SONLIB} &&\
    cd sonLib &&\
    git checkout ${VERSION_SONLIB} &&\
    make &&\

### Installing crossmap
    cd /opt &&\
    git clone ${URL_HAL} &&\
    cd hal &&\
    git checkout ${VERSION_HAL_COMMIT} &&\
    make &&\
    strip /opt/hal/bin/*; true



