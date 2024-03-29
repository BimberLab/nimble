FROM centos:7

# Update mirrors
RUN yum update -y

# Set proper encoding for pip 20.2+
ENV LANG en_US.utf8
ENV LC_ALL en_US.utf8

# Install python3, pip, and dependencies
RUN yum group install -y "Development Tools" && \
    yum install -y ncurses-devel bzip2-devel xz-devel zlib-devel wget glibc-devel python3-devel && \
    yum install -y python3 python-pip && \
    yum clean all && \
    rm -rf /var/cache/yum

# Install samtools
RUN cd /tmp &&\
        wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 &&\
        tar xvjf samtools-1.13.tar.bz2 &&\
        cd samtools-1.13 &&\
        ./configure --prefix=/usr/local &&\
        make &&\
        make install

# Install nimble
ADD . /nimble
RUN ["pip3", "install", "./nimble"]

# Download the latest aligner version
RUN ["python3", "-m", "nimble", "download"]

# Install jemalloc with profiling from the stable-4 branch
RUN cd /opt \
    && git clone --branch stable-4 https://github.com/jemalloc/jemalloc.git \
    && cd jemalloc \
    && ./autogen.sh \
    && ./configure --prefix=/usr/local --enable-prof \
    && make \
    && make install_bin install_include install_lib

ENV LD_PRELOAD="/usr/local/lib/libjemalloc.so"
#ENV MALLOC_CONF=prof_leak:true,lg_prof_sample:19,prof_final:true,prof_prefix:/work/jemalloc_profile

ENTRYPOINT ["python3", "-m", "nimble"]
