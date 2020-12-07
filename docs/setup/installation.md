# Installation 

You can install RAPIDS using Docker (the fastest), or native instructions for MacOS and Ubuntu

=== "Docker"
    
    1.  Install [Docker](https://docs.docker.com/desktop/)

    2.  Pull our RAPIDS container
        ``` bash
        docker pull agamk/rapids:latest`
        ```

    3.  Run RAPIDS\' container (after this step is done you should see a
        prompt in the main RAPIDS folder with its python environment active)

        ``` bash
        docker run -it agamk/rapids:latest
        ```

    4.  Pull the latest version of RAPIDS

        ``` bash
        git pull
        ```
    
    5. Make RAPIDS script executable
        ```bash
        chmod +x rapids
        ```
    
    6.  Check that RAPIDS is working
        ``` bash
        ./rapids -j1
        ```
    7.  *Optional*. You can edit RAPIDS files with `vim` but we recommend using `Visual Studio Code` and its `Remote Containers` extension

        ??? info "How to configure Remote Containers extension"

            - Make sure RAPIDS container is running
            - Install the [Remote - Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
            - Go to the `Remote Explorer` panel on the left hand sidebar
            - On the top right dropdown menu choose `Containers`
            - Double click on the `agamk/rapids` container in the`CONTAINERS` tree
            - A new VS Code session should open on RAPIDS main folder insidethe container.

=== "MacOS"
    We tested these instructions in Catalina

    1.  Install [brew](https://brew.sh/)

    2.  Install MySQL

        ``` bash
        brew install mysql
        brew services start mysql
        ```

    3.  Install R 4.0, pandoc and rmarkdown. If you have other instances of R, we recommend uninstalling them

        ``` bash
        brew install r
        brew install pandoc
        Rscript --vanilla -e 'install.packages("rmarkdown", repos="http://cran.us.r-project.org")'
        ```

    4.  Install miniconda (restart your terminal afterwards)

        ``` bash
        brew cask install miniconda
        conda init zsh # (or conda init bash)
        ```

    5.  Clone our repo

        ``` bash
        git clone https://github.com/carissalow/rapids
        ```

    6.  Create a python virtual environment

        ``` bash
        cd rapids
        conda env create -f environment.yml -n rapids
        conda activate rapids
        ```

    7.  Install R packages and virtual environment:

        ``` bash
        snakemake -j1 renv_install
        snakemake -j1 renv_restore
           
        ```

        !!! note
            This step could take several minutes to complete, especially if you have less than 3Gb of RAM or packages need to be compiled from source. Please be patient and let it run until completion.
    
    5. Make RAPIDS script executable
        ```bash
        chmod +x rapids
        ```

    8.  Check that RAPIDS is working
        ``` bash
        ./rapids -j1
        ```

=== "Ubuntu"

    We tested on Ubuntu 18.04 & 20.04

    1.  Install dependencies

        ``` bash
        sudo apt install libcurl4-openssl-dev
        sudo apt install libssl-dev
        sudo apt install libxml2-dev
        ```

    2.  Install MySQL

        ``` bash
        sudo apt install libmysqlclient-dev
        sudo apt install mysql-server
        ```

    3.  Add key for R's repository.

        ``` bash
        sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
        ```

    4. Add R's repository

        1. For 18.04    
        ``` bash
        sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
        ```

        1. For 20.04
        ``` bash
        sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
        ```

    5. Install R 4.0. If you have other instances of R, we recommend uninstalling them

        ``` bash
        sudo apt update
        sudo apt install r-base
        ```

    6.  Install Pandoc and rmarkdown

        ``` bash
        sudo apt install pandoc
        Rscript --vanilla -e 'install.packages("rmarkdown", repos="http://cran.us.r-project.org")'
        ```

    7.  Install git

        ``` bash
        sudo apt install git
        ```

    8.  Install [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

    9.  Restart your current shell

    10. Clone our repo:

        ``` bash
        git clone https://github.com/carissalow/rapids
        ```

    11. Create a python virtual environment:

        ``` bash
        cd rapids
        conda env create -f environment.yml -n MY_ENV_NAME
        conda activate MY_ENV_NAME
        ```

    7.  Install R packages and virtual environment:

        ``` bash
        snakemake -j1 renv_install
        snakemake -j1 renv_restore
           
        ```

        !!! note
            This step could take several minutes to complete, especially if you have less than 3Gb of RAM or packages need to be compiled from source. Please be patient and let it run until completion.

    5. Make RAPIDS script executable
        ```bash
        chmod +x rapids
        ```

    8.  Check that RAPIDS is working
        ``` bash
        ./rapids -j1
        ```
