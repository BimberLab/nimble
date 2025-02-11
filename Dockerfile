FROM python:3.9-slim

RUN apt-get update && apt-get install -y \
    samtools \
    wget \
    git \
    autoconf \
    automake \
    gcc \
    make \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --branch stable-4 https://github.com/jemalloc/jemalloc.git /opt/jemalloc && \
    cd /opt/jemalloc && \
    ./autogen.sh && \
    ./configure --prefix=/usr/local --enable-prof && \
    make && \
    make install_bin install_include install_lib && \
    rm -rf /opt/jemalloc

RUN apt-get remove -y autoconf automake gcc make git wget && apt-get autoremove -y

ENV LD_PRELOAD="/usr/local/lib/libjemalloc.so"

ADD . /nimble
RUN pip install --no-cache-dir ./nimble

RUN python3 -m nimble download

ENTRYPOINT ["python3", "-m", "nimble"]