#!/usr/bin/env bash

UBUNTU_VERSION=${1:-"18.04"}
PYTHON_VERSION=${2:-"3.6.0"}
CWLTOOL_VERSION=${3:-"main"}
CWLTEST_VERSION=${4:-"main"}

WORKING_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
echo "Running conformance tests with Python ${PYTHON_VERSION} in dockerized Ubuntu $UBUNTU_VERSION"
echo "Using cwltool==${CWLTOOL_VERSION} and cwltest==${CWLTEST_VERSION}"
echo "Working directory $WORKING_DIR"
docker rmi cwltool:latest --force
docker build --no-cache --build-arg UBUNTU_VERSION=$UBUNTU_VERSION \
                        --build-arg PYTHON_VERSION=$PYTHON_VERSION \
                        --build-arg CWLTOOL_VERSION=$CWLTOOL_VERSION \
                        --build-arg CWLTEST_VERSION=$CWLTEST_VERSION \
                        --rm -t cwltool:latest .
docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock \
                    -v ${WORKING_DIR}:${WORKING_DIR} \
                    --workdir ${WORKING_DIR} \
                    cwltool:latest \
                    ${WORKING_DIR}/run_test.sh --junit-xml=result.xml -j2 RUNNER=cwltool \
                    EXTRA="--tmpdir-prefix=${WORKING_DIR}/temp/tmp --tmp-outdir-prefix=${WORKING_DIR}/temp/tmp"