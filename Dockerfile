FROM continuumio/miniconda:4.5.12

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
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && \
    echo "source activate scsim" > ~/.bashrc 

# Install monovar conda environment
COPY monovar-env.yml /tmp/monovar-env.yml
RUN conda env create -f /tmp/monovar-env.yml

# Copy scsim into container
COPY src /opt/scsim

ENV PATH="/opt/conda/envs/scsim/bin:/opt/dwgsim:${PATH}"


ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
