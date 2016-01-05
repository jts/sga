# SGA docker
# Version 0.1
# Runs SGA ( http://github.com/jts/sga ) 
#
# Creates SGA-in-a-docker, such that you can run the dockerized SGA
# much as you would run the local SGA executable.
# 
# To build:
#    docker build -t sga-docker [path containing DOCKERFILE]
# if the Dockerfile is local, otherwise:
#    docker build -t sga-docker github.com/jts/sga
#    
# Once built, you can invoke sga-docker as you would sga, eg:
#    docker run -v /path/to/local/data/data/:/data/ \
#       sga-docker somatic-variant-filters --annotate-only --threads=2 \
#        --tumor-bam=/data/tumor.bam --normal-bam=/data/normal.bam \
#        --reference=/data/reference.fa /data/input.vcf
#
# You can also invoke other sga scripts, such as sga-preqc-report:
#   docker run -v /path/of/sga/repo:/data/ \
#       --entrypoint /src/sga-0.10.14/src/bin/sga-preqc-report.py \
#       sga-docker /data/src/examples/preqc/{bird,fish,human,oyster,snake,yeast}.preqc
#
FROM ubuntu:14.04
MAINTAINER Jonathan Dursi <Jonathan.Dursi@oicr.on.ca> 
LABEL Description="Runs SGA variant annotation on candidate indels against tumour and normal bams" Vendor="OICR" Version="0.1"
VOLUME /data
WORKDIR /data

# get ubuntu packages
# incl python matplotlib for preqc
RUN apt-get update && \
    apt-get install -y \
        automake \
        autotools-dev \
        build-essential \
        cmake \
        libhts-dev \
        libhts0 \
        libjemalloc-dev \
        libsparsehash-dev \
        libz-dev \
        python-matplotlib \
        wget \
        zlib1g-dev 

# build remaining dependencies:
# bamtools
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/pezmaster31/bamtools/archive/v2.4.0.tar.gz && \
    tar -xzvf v2.4.0.tar.gz && \
    rm v2.4.0.tar.gz && \
    cd bamtools-2.4.0 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# build SGA
RUN mkdir -p /src && \
    cd /src && \
    wget https://github.com/jts/sga/archive/v0.10.14.tar.gz && \
    tar -xzvf v0.10.14.tar.gz && \
    rm v0.10.14.tar.gz && \
    cd sga-0.10.14/src && \
    ./autogen.sh && \
    ./configure --with-bamtools=/deps/bamtools-2.4.0 --with-jemalloc=/usr --prefix=/usr/local && \
    make && \
    make install


ENTRYPOINT ["/usr/local/bin/sga"]
CMD ["--help"]
