.. _install-page:

Installation
===============

This instructions have been tested on macOS (Catalina and Mojave) and Ubuntu 16.04. If you find a problem, please report it on our GitHub page.

macOS (tested on Catalina 10.15)
--------------------------------

#. Install dependenies (Homebrew if not installed):

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
    - ``conda env create -f environment.yml -n rapids``
    - ``conda activate rapids``

#. Install r packages and virtual environment:

    - ``snakemake packrat_install``
    - ``snakemake packrat_init``
    - ``snakemake packrat_restore``
        - This step will take several minutes to complete. Please be patient and let it run until completion. 

#. See Usage section below. 


Linux (tested on Ubuntu 16.04)
------------------------------

#. Install dependenies (Homebrew - if not installed):

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
        - This step will take several minutes to complete. Please be patient and let it run until completion. 

#. See Usage section below.


Usage
======
Once you have the installation for your specific operating system complete, you can follow these steps to get starting using the pipeline.

#. Configure the participants you want to analyze:

    - Create one file per participant in the ``rapids/data/external/`` directory. The file should NOT have an extension (i.e. no .txt). The name of the file will become the label for that participant in the pipeline. 
    - The first line of the file should be a comma separated list with each of the device ID numbers for that participant as it appears in AWARE.
        - If AWARE is removed and reinstalled on the device, a new device ID is generated.
    - The second line should list the device's operating system (Android or iOS)
    - As an example. Let's say participant `p00` had 2 AWARE device_id numbers and was running Android OS. Their file would be named `p00` and contain:

        .. code-block:: bash

            3a7b0d0a-a9ce-4059-ab98-93a7b189da8a,44f20139-50cc-4b13-bdde-0d5a3889e8f9
            android


#. Configure the database connection:

    - Create an empty file called `.env` in the root directory (rapids/)
    - Add and complete the following lines:

        | ``[MY_GROUP_NAME]``
        | ``user=MyUSER``
        | ``password=MyPassword``
        | ``host=MyIP``
        | ``port=3306``
        - Replace your database specific credentials with those listed above.  


#. Once the all of the installation and configurations has been completed the following command can be run to pull the default sample dataset that comes with this project.::

    $ snakemake


This pulls sample data from AWARE_ and processes it with the default rules that come with RAPIDS.

.. _the-install-note:

.. note::
    - Ensure that ``MY_GROUP_NAME`` is the same value for GROUP in the ``DATABASE_GROUP`` variable in the ``config.yaml`` file located in the root directory (rapids/config.yaml). 
    - Ensure that your list of ``SENSORS`` in the ``config.yaml`` file correspond to the sensors used in the ``all`` rule in the ``Snakefile`` file (See :ref:`rapids-structure` for more information)
    



.. _bug: https://github.com/Homebrew/linuxbrew-core/issues/17812
.. _instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
.. _brew: https://docs.brew.sh/Homebrew-on-Linux
.. _AWARE: https://awareframework.com/what-is-aware/