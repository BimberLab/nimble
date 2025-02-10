FROM centos:7

# Set the mirrors to vault, as centos7 is EoL
RUN sed -i s/mirror.centos.org/vault.centos.org/g /etc/yum.repos.d/*.repo && \
    sed -i s/^#.*baseurl=http/baseurl=http/g /etc/yum.repos.d/*.repo && \
    sed -i s/^mirrorlist=http/#mirrorlist=http/g /etc/yum.repos.d/*.repo

# Update mirrors
RUN yum update -y

# Set proper encoding for pip 20.2+
ENV LANG=en_US.utf8
ENV LC_ALL=en_US.utf8

# Install developer tools and dependencies
RUN yum group install -y "Development Tools" && \
    yum install -y ncurses-devel bzip2-devel xz-devel zlib-devel wget glibc-devel python3-devel \
                   openssl11-devel epel-release.noarch libffi-devel centos-release-scl devtoolset-8

# Install OpenSSL 1.1.1k
RUN wget https://openssl.org/source/openssl-1.1.1k.tar.gz && \
    tar -xzvf openssl-1.1.1k.tar.gz && \
    cd openssl-1.1.1k && \
    ./config --prefix=/usr --openssldir=/etc/ssl --libdir=lib no-shared zlib-dynamic && \
    make && \
    make install

# Install Python 3.9.6 with OpenSSL 1.1.1 support
RUN wget https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz && \
    tar xzf Python-3.9.6.tgz && \
    cd Python-3.9.6 && \
    sed -i 's/PKG_CONFIG openssl /PKG_CONFIG openssl11 /g' configure && \
    ./configure --enable-optimizations && \
    make install && \
    pip3 install --upgrade pip

# Clean up
RUN yum clean all && rm -rf /var/cache/yum

# Install samtools and create tmp folders:
RUN cd /tmp &&\
        wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 &&\
        tar xvjf samtools-1.21.tar.bz2 &&\
        cd samtools-1.21 &&\
        ./configure --prefix=/usr/local &&\
        make &&\
        make install &&\
        # This is to avoid the numba 'cannot cache function' error, such as: https://github.com/numba/numba/issues/5566
        mkdir /numba_cache && chmod -R 777 /numba_cache &&\
        mkdir /mpl_cache && chmod -R 777 /mpl_cache
        
# NOTE: this is required when running as non-root. Setting MPLCONFIGDIR removes a similar warning.
ENV NUMBA_CACHE_DIR=/numba_cache
ENV MPLCONFIGDIR=/mpl_cache
ENV CFLAGS="-std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE"
ENV LDFLAGS="-lrt"

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
