FROM continuumio/miniconda:4.6.14

RUN apt-get update --fix-missing && \
    apt-get install -y \
        build-essential \
        zlib1g-dev \
        libncurses-dev \
        vim \
        wget


# Install dwgsim
RUN git clone --recursive https://github.com/nh13/DWGSIM.git /opt/dwgsim && \
    cd /opt/dwgsim && \
    git submodule init && \
    git submodule update && \
    make --directory /opt/dwgsim

# Install scsim conda environment and make default
COPY envs/scsim-env.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
#RUN ln -s /opt/conda//etc/profile.d/conda.sh /etc/profile.d/conda.sh


# Install monovar conda environment
COPY envs/monovar-env.yml /tmp/monovar-env.yml
#RUN conda env create -f /tmp/monovar-env.yml

# Install bam-readcount conda environment
COPY envs/bam-readcount-env.yml /tmp/bam-readcount-env.yml
#RUN conda env create -f /tmp/bam-readcount-env.yml

# Copy scsim into container
COPY src /opt/scsim/src

ENV PATH="/opt/conda/envs/scsim/bin:/opt/dwgsim:${PATH}"

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
