FROM continuumio/miniconda:4.6.14

RUN apt-get update --fix-missing && \
    apt-get install -y \
        build-essential \
        zlib1g-dev \
        libncurses-dev \
        vim \
        wget \
        bash


# Install dwgsim
RUN git clone --recursive https://github.com/nh13/DWGSIM.git /opt/dwgsim && \
    cd /opt/dwgsim && \
    git submodule init && \
    git submodule update && \
    make --directory /opt/dwgsim

# Install scsim conda environment and make default
COPY envs/scsim-env.yml /envs/environment.yml
RUN conda env create -f /envs/environment.yml
RUN rm /etc/profile.d/conda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh


# Install monovar conda environment
COPY envs/monovar-env.yml /envs/monovar-env.yml
#RUN conda env create -f /envs/monovar-env.yml

# Install bam-readcount conda environment
COPY envs/bam-readcount-env.yml /envs/bam-readcount-env.yml
#RUN conda env create -f /envs/bam-readcount-env.yml

# Copy scsim into container
COPY src /src

ENV PATH="/opt/conda/envs/scsim/bin:/opt/dwgsim:${PATH}"

RUN mkdir /.conda && chmod -R a+w /.conda

ENTRYPOINT [ "/usr/bin/tini", "--" ]
