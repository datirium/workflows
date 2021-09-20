
#####################################################################################
# Docker image for running conformance tests                                        #
# with custom UBUNTU_VERSION, PYTHON_VERSION, CWLTOOL_VERSION, and  CWLTEST_VERSION #
#####################################################################################
# Build Cmd:        docker build --no-cache --rm -t cwltool:latest .                #
#####################################################################################

# can be provided through --build-arg PARAM=value
ARG UBUNTU_VERSION="20.04"
ARG PYTHON_VERSION="3.8.10"
ARG CWLTOOL_VERSION="main"
ARG CWLTEST_VERSION="main"

FROM ubuntu:${UBUNTU_VERSION}
LABEL maintainer="heylel.b.sh@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /tmp

ARG PYTHON_VERSION
ARG CWLTOOL_VERSION
ARG CWLTEST_VERSION

ENV PYTHON_URL "https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz"
ENV CWLTOOL_URL "https://github.com/common-workflow-language/cwltool.git"
ENV CWLTEST_URL "https://github.com/common-workflow-language/cwltest.git"

RUN echo "Installing dependencies" && \
    apt-get update && \
    apt-get install -y gcc build-essential \
                       git wget curl zlib1g-dev libmysqlclient-dev libffi-dev libssl-dev \
                       ca-certificates \
                       nodejs mysql-client apt-transport-https libsqlite3-dev \
                       gnupg-agent software-properties-common && \
    echo "Installing Python" && \
    wget ${PYTHON_URL} && \
    tar xzf Python-${PYTHON_VERSION}.tgz && \
    cd Python-${PYTHON_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    pip3 install -U pip && \
    echo "Installing docker-ce-cli" && \
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add - && \
    add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" && \
    apt-get update && \
    apt-get -y install docker-ce-cli && \
    echo "Installing cwltool" && \
    git clone ${CWLTOOL_URL} && \
    cd cwltool && \
    git checkout ${CWLTOOL_VERSION} && \
    pip3 install . && \
    cd .. && \
    echo "Installing cwltest" && \
    git clone ${CWLTEST_URL} && \
    cd cwltest && \
    git checkout ${CWLTEST_VERSION} && \
    pip3 install . && \
    cd .. && \
    # cleaning up
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true