FROM ubuntu:21.10

# See: https://stackoverflow.com/questions/44331836/apt-get-install-tzdata-noninteractive
ENV DEBIAN_FRONTEND=noninteractive

#samtools' dependencies
RUN apt-get update -y \
		&& apt-get install -y libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev python3-pip wget \
		&& apt-get clean

#samtools
RUN cd /tmp \
		&& wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 \
		&& tar xvjf samtools-1.13.tar.bz2 \
		&& rm samtools-1.13.tar.bz2 \
		&& cd samtools-1.13 \
		&& ./configure --prefix=/usr/local \
		&& make \
		&& make install

# Install nimble
ADD . /nimble
RUN ["pip", "install", "./nimble"]

# Download the latest aligner version
RUN ["python3", "-m", "nimble", "download"]

ENTRYPOINT ["python3", "-m", "nimble"]
