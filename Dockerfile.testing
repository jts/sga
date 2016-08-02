#
# Dockerfile meant for build testing only.
#
# Synopsis:
# - install prerequisites from debian:stable
# - install gcc-6 from debian:testing
# - clone bamtools (NOTE: fix for gcc-6 not included in currently latest tag 2.4.0)
# - build bamtools
# - install bamtools in /usr/local
# - copy (sga/) src from Dockerfile context
# - build sga
# - install sga in /usr/local
#
FROM debian:stable
MAINTAINER Matei David <matei.david.at.oicr.on.ca>
ARG DEBIAN_FRONTEND=noninteractive
LABEL Description="Dockerfile meant for build testing only."
WORKDIR /tmp

# enable debian:testing
RUN echo 'APT::Default-Release "stable";' >/etc/apt/apt.conf.d/99defaultrelease && \
    mv /etc/apt/sources.list /etc/apt/sources.list.d/stable.list && \
    sed 's/stable/testing/g' </etc/apt/sources.list.d/stable.list >/etc/apt/sources.list.d/testing.list

# install prerequisites from stable
RUN apt-get update -y && \
    apt-get install -y \
        automake \
        autotools-dev \
        build-essential \
        cmake \
        git \
        libjemalloc-dev \
        libsparsehash-dev \
        libz-dev \
        wget

# install gcc-6 from testing
RUN apt-get install -y -t testing \
    g++-6
ENV CC=gcc-6
ENV CXX=g++-6

# build bamtools
RUN git clone --depth 1 https://github.com/pezmaster31/bamtools.git && \
    cd bamtools && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

ADD src /tmp/sga/src/
RUN cd /tmp/sga/src && \
    ./autogen.sh && \
    ./configure --with-bamtools=/usr/local --with-jemalloc=/usr && \
    make && \
    make install

VOLUME /data
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/sga"]
CMD ["--help"]
