# Installation 

You can install RAPIDS using Docker (the fastest), or native instructions for MacOS and Linux (Ubuntu). Windows is supported through Docker or WSL.

=== "Docker"
    
    1.  Install [Docker](https://docs.docker.com/desktop/)

    2.  Pull our RAPIDS container
        ``` bash
        docker pull moshiresearch/rapids:latest
        ```

    3.  Run RAPIDS\' container (after this step is done you should see a
        prompt in the main RAPIDS folder with its python environment active)

        ``` bash
        docker run -it moshiresearch/rapids:latest
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

        ??? info "How to configure the Remote Containers extension"
            
            - Make sure RAPIDS Docker container is running

            - Install VS Code and its [Remote - Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

            - Click the `Remote Explorer` icon on the left-hand sidebar (the icon is a computer monitor)

            - On the top right dropdown menu, choose `Containers`

            - Right-click on the `moshiresearch/rapids` container in the `CONTAINERS` tree and select `Attach to Container`. A new VS Code window should open

            - In the new window, open the `/rapids/` folder via the `File/Open...` menu

            - Run RAPIDS inside a terminal in VS Code. Open one with the `Terminal/New Terminal` menu

    !!! warning
        If you installed RAPIDS using Docker for Windows on Windows 10, the container will have [limits](https://stackoverflow.com/questions/43460770/docker-windows-container-memory-limit) on the amount of RAM it can use. If you find that RAPIDS crashes due to running out of memory, [increase](https://stackoverflow.com/a/56583203/6030343) this limit.

=== "MacOS"
    We tested these instructions in Catalina and Big Sur

    ??? info "M1 Macs"
        RAPIDS can run on M1 Macs, the only changes as of Mar 17, 2022 are:

        - Brew and everything installed with it needs to be setup under Rosetta (x86 arch) due to incompatibility issues with some R libraries and python packages. To do this, run your terminal [via Rosetta](https://www.youtube.com/watch?v=nv2ylxro7rM&t=138s), then proceed with our installation commands. 
        - There is a bug related to timezone codes. We set the correct `TZ_DIR` in `renv/activate.R` (line #19) `Sys.setenv("TZDIR" = file.path(R.home(), "share", "zoneinfo"))` (RAPIDS does this automatically).

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
        brew install --cask miniconda
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

    We tested RAPIDS on Ubuntu 18.04 & 20.04. Note that the necessary Python and R packages are available in other Linux distributions, so if you decide to give it a try, let us know and we can update these docs.

    1.  Install dependencies

        ``` bash
        sudo apt install libcurl4-openssl-dev
        sudo apt install libssl-dev
        sudo apt install libxml2-dev
        sudo apt install libglpk40
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

        === "Ubuntu 18.04 Bionic"

            ``` bash
            sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
            ```

        === "Ubuntu 20.04 Focal"

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

    7.  Install the R virtual environment management package (renv)

        ``` bash
        snakemake -j1 renv_install
        ```

    8. Restore the R virtual environment

        === "Ubuntu 18.04 Bionic (fast)"

            Run the following command to restore the R virtual environment using [RSPM](https://packagemanager.rstudio.com/client/#/repos/1/overview) binaries
            ```bash
            R -e 'renv::restore(repos = c(CRAN = "https://packagemanager.rstudio.com/all/__linux__/bionic/latest"))'
            ```

        === "Ubuntu 20.04 Focal (fast)"

            Run the following command to restore the R virtual environment using [RSPM](https://packagemanager.rstudio.com/client/#/repos/1/overview) binaries
            ```bash
            R -e 'renv::restore(repos = c(CRAN = "https://packagemanager.rstudio.com/all/__linux__/focal/latest"))'
            ```
            
        === "Ubuntu (slow)"

            If the fast installation command failed for some reason, you can restore the R virtual environment from source:
            ```bash
            R -e 'renv::restore()'
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
    
=== "Windows"

    There are several options varying in complexity:

    - You can use our Docker instructions (tested)
    - You can use our Ubuntu 20.04 instructions on [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install-win10) (not tested but it will likely work)
    - Native installation (experimental). If you would like to contribute to RAPIDS you could try to install MySQL, miniconda, Python, and R 4.0+ in Windows and restore the Python and R virtual environments using steps 6 and 7 of the instructions for Mac. You can [get in touch](../../team) if you would like to discuss this with the team.
