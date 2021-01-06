# getting base image ubuntu
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libmysqlclient-dev \
    libglpk40 \
    mysql-server
RUN apt-get update && apt-get install -y gnupg
RUN apt-get update && apt-get install -y software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt update && apt install -y r-base
RUN apt install -y pandoc
RUN apt install -y git
RUN apt-get update && apt-get install -y vim
RUN apt-get update && apt-get install -y nano
RUN apt update && apt install -y unzip
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
RUN git clone https://github.com/carissalow/rapids
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
RUN conda update -n base -c defaults conda
WORKDIR /rapids
RUN conda env create -f environment.yml -n rapids
RUN Rscript --vanilla -e 'install.packages("rmarkdown", repos="http://cran.us.r-project.org")'
RUN R -e 'renv::restore(repos = c(CRAN = "https://packagemanager.rstudio.com/all/__linux__/focal/latest"))'
ADD https://osf.io/587wc/download data/external
RUN mv data/external/download data/external/rapids_example.sql.zip
RUN unzip data/external/rapids_example.sql.zip
RUN cp rapids_example.sql data/external/rapids_example.sql
RUN rm data/external/rapids_example.sql.zip
RUN rm rapids_example.sql
RUN echo "source activate rapids" > ~/.bashrc
ENV PATH /opt/conda/envs/rapids/bin:$PATH