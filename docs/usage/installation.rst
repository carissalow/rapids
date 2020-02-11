Installation
===============

This instructions have been tested on MacOS Catalina and Mojave and Ubuntu 16.04. If you find a problem, please report it.

Mac OS (tested on Catalina)
----------------------------

#. Install dependenies (if not installed):

    - Install brew_ for Mac: ``/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"``

#. Install MySQL

    - ``brew install mysql``
    - ``brew services start mysql``

#. Install R, pandoc and rmarkdown:

    - ``brew install r``
    - ``brew install pandoc``
    - ``R -e 'install.packages(c( "rmarkdown"), repos = "http://cran.us.r-project.org")'``

#. Install miniconda:

    - ``brew cask install miniconda``
    - ``conda init zsh`` or ``conda init bash``
    - Restart terminal if necessary

#. Clone our repo:

    - ``git clone https://github.com/carissalow/rapids``

#. Create a python virtual environment:

    - ``cd rapids``
    - ``conda env create -f environment.yml -n MY_ENV_NAME``
    - ``conda activate MY_ENV_NAME``

#. Install r packages and virtual environment:

    - ``snakemake packrat_install``
    - ``snakemake packrat_init``
    - ``snakemake packrat_restore``

#. Configure the participants to analyze:

    - Create a file per participant in the ``/data/external`` folder, no extension is necessary, its name will be the label for that participant in the pipeline: ``/data/external/pxx``
    - Add a line with the device_id(s) of that participant as it appears on the database. If multiple device ids, all data for this participant will be relabeled with the last one
    - Add a line with the mobile platform (android, or ios)
    - For example:

        .. code-block:: bash

            3a7b0d0a-a9ce-4059-ab98-93a7b189da8a,44f20139-50cc-4b13-bdde-0d5a3889e8f9
            android

#. Configure the db connection:

    - Create an empty .env file in the root folder
    - Add and complete the following lines:

        | ``[MY_GROUP_NAME]``
        | ``user=MyUSER``
        | ``password=MyPassword``
        | ``host=MyIP``
        | ``port=3306``

Linux (tested on Ubuntu 16.04)
------------------------------

#. Install dependenies (if not installed):

    - ``sudo apt-get install libmariadb-client-lgpl-dev libxml2-dev libssl-dev``
    - Install brew_ for linux and add the following line to ~/.bashrc: ``export PATH=$HOME/.linuxbrew/bin:$PATH``
    - ``source ~/.bashrc``

#. Install MySQL

    - ``brew install mysql``
    - ``brew services start mysql``

#. Install R, pandoc and rmarkdown:

    - ``brew install r``
    - ``brew install gcc@6`` (needed due to this bug_)
    - ``HOMEBREW_CC=gcc-6 brew install pandoc``
    - ``R -e 'install.packages(c( "rmarkdown"), repos = "http://cran.us.r-project.org")'``

#. Install miniconda using these instructions_

#. Clone our repo:

    - ``git clone https://github.com/carissalow/rapids``

#. Create a python virtual environment:

    - ``cd rapids``
    - ``conda env create -f environment.yml -n MY_ENV_NAME``
    - ``conda activate MY_ENV_NAME``

#. Install r packages and virtual environment:

    - ``snakemake packrat_install``
    - ``snakemake packrat_init``
    - ``snakemake packrat_restore``

#. Configure the participants to analyze:

    - Create a file per participant in the ``/data/external`` folder, no extension is necessary, its name will be the label for that participant in the pipeline: ``/data/external/pxx``
    - Add a line with the device_id(s) of that participant as it appears on the database. If multiple device ids, all data for this participant will be relabeled with the last one
    - Add a line with the mobile platform (android, or ios)
    - For example:

        .. code-block:: bash

            3a7b0d0a-a9ce-4059-ab98-93a7b189da8a,44f20139-50cc-4b13-bdde-0d5a3889e8f9
            android

#. Configure the db connection:

    - Create an empty .env file in the root folder
    - Add and complete the following lines:

        | ``[MY_GROUP_NAME]``
        | ``user=MyUSER``
        | ``password=MyPassword``
        | ``host=MyIP``
        | ``port=3306``


.. note::
    - Ensure that ``MY_GROUP_NAME`` is the same value for GROUP in the ``DATABASE_GROUP`` variable in the config.yaml file. 
    - Ensure that your list of ``SENSORS`` in the config.yaml file correspond to the sensors used in all rules in the Snakefile file (See Pipeline Structure Section for more information TBD)



.. _bug: https://github.com/Homebrew/linuxbrew-core/issues/17812
.. _instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
.. _brew: https://docs.brew.sh/Homebrew-on-Linux