FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update


# -----------------------------------------

RUN apt-get -y update \
	&& apt-get -y install gcc \
	&& apt-get -y install g++ \
	&& apt-get -y install autoconf \
	&& apt-get -y install zlib1g-dev \
	&& apt-get -y install wget \
	&& apt-get -y install pkg-config

RUN pip install --upgrade pip \
    && pip install -q pyvcf

RUN git clone https://github.com/vcftools/vcftools.git \
    && cd vcftools \
    && ./autogen.sh \
    && ./configure \
    && make \
    && make install

RUN wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.1/vcf_validator_linux \
    && chmod 755 vcf_validator_linux \
    && mv vcf_validator_linux /kb/deployment/bin 

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
