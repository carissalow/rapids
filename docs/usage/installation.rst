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

#. Install R packages and virtual environment:

    - ``snakemake renv_install``
    - ``snakemake renv_init``
    - ``snakemake renv_restore``
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

#. Install R packages and virtual environment:

    - ``snakemake renv_install``
    - ``snakemake renv_init``
    - ``snakemake renv_restore``
        - This step will take several minutes to complete. Please be patient and let it run until completion. 

#. See Usage section below.


Usage
======
Once you have the installation for your specific operating system complete, you can follow these steps to start using RAPIDS.

.. _db-configuration:

#. Configure the database connection:

    - Create an empty file called `.env` in the root directory (``rapids/``)
    - Add the following lines and replace your database specific credentials (user, password, and host):

        .. code-block:: bash
        
            [MY_GROUP]
            user=MyUSER
            password=MyPassword
            host=MyIP
            port=3306
            database=MyDB

        .. note::

            ``MY_GROUP`` is a custom label you assign when setting up the database configuration. It has to match ``DATABASE_GROUP`` in the ``config.yaml`` file_. It does not have to relate to your database credentials.

#. Configure the participants you want to analyze:

    - **Automatically**. You can automatically include all devices that are stored in the ``aware_device`` table, if you have especial requirements see the Manual configuration::

        snakemake download_participants

    - **Manually**. Create one file per participant in the ``rapids/data/external/`` directory. The file should NOT have an extension (i.e. no .txt). The name of the file will become the label for that participant in the pipeline.

        - The first line of the file should be the Aware ``device_id`` for that participant. If one participant has multiple device_ids (i.e. Aware had to be re-installed), add all device_ids separated by commas.
        - The second line should list the device's operating system (``android`` or ``ios``)
        - The third line is a human friendly label that will appear in any plots for that participant.
        - The forth line contains a start and end date separated by a comma (e.g. ``20201301,20202505``). Only data wihtin these dates will be included in the pipeline.

    For example, let's say participant `p01` had two AWARE device_ids and they were running Android between Feb 1st 2020 and March 3rd 2020. Their participant file would be named ``p01`` and contain:

        .. code-block:: bash

            3a7b0d0a-a9ce-4059-ab98-93a7b189da8a,44f20139-50cc-4b13-bdde-0d5a3889e8f9
            android
            Participant01
            2020/02/01,2020/03/03

#. Configure the sensors to process:

    - The variable ``SENSORS`` in the ``config.yaml`` file_ should match existent sensor tables in your Aware database (See :ref:`rapids-structure` for more information). Each item in this list will be processed in RAPIDS.

    .. note::

        It is beneficial to list all collected sensors even if you don't plan to include them in a model later on in the pipeline. This is because we use all data available to estimate whether the phone was sensing data or not (i.e. to know if Aware crashed or the battery died). See :ref:`PHONE_VALID_SENSED_DAYS<phone-valid-sensed-days>` for more information.

#. Execute RAPIDS

    - Standard execution over a single core::

        snakemake -j1
    
    - Standard execution over multiple cores::

        snakemake -j8

    - Force a rule (useful if you modify your code and want to update its results)::

        snakemake -j1 -R RULE_NAME

.. _bug: https://github.com/Homebrew/linuxbrew-core/issues/17812
.. _instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
.. _brew: https://docs.brew.sh/Homebrew-on-Linux
.. _AWARE: https://awareframework.com/what-is-aware/
.. _file: https://github.com/carissalow/rapids/blob/master/config.yaml#L22
