.. _install-page:

Installation
===============

These instructions have been tested on macOS (Catalina and Mojave) and Ubuntu 16.04. If you find a problem, please create a GitHub issue or contact us. If you want to test RAPIDS quickly try our docker image or follow the Linux instructions on a virtual machine.

Docker (the fastest and easiest way)
------------------------------------

#. Install docker

#. Pull RAPIDS' container

    ``docker pull agamk/rapids:latest``

#. Run RAPIDS' container (after this step is done you should see a prompt in the main RAPIDS folder with its python environment active)

    ``docker run -it agamk/rapids:latest``

#. Optional. You can start editing files with vim but we recommend using Visual Studio Code and its Remote extension

    - Make sure RAPIDS container is running
    - Install the Remote - Containers extension: https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers
    - Go to the ``Remote Explorer`` panel on the left hand sidebar
    - On the top right dropdown menu choose ``Containers``
    - Double click on the ``agamk/rapids`` container in the ``CONTAINERS`` tree
    - A new VS Code session should open on RAPIDS main folder inside the container.

#. See Usage section below.


macOS (tested on Catalina 10.15)
--------------------------------

#. Install dependencies (Homebrew if not installed):

    - Install brew_ for Mac: ``/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"``

#. Install MySQL

    - ``brew install mysql``
    - ``brew services start mysql``

#. Install R 4.0 and pandoc. If you have other instances of R, we recommend uninstalling them.

    - ``brew install r``
    - ``brew install pandoc``

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

    - ``snakemake -j1 renv_install``
    - ``snakemake -j1 renv_restore``

      - This step could take several minutes to complete, especially if you have less than 3Gb of RAM or packages need to be compiled from source. Please be patient and let it run until completion.  

#. See Usage section below. 


Linux (tested on Ubuntu 18.04 & 20.04)
---------------------------------------

#. Install dependencies :

    - ``sudo apt install libcurl4-openssl-dev``
    - ``sudo apt install libssl-dev``
    - ``sudo apt install libxml2-dev``

#. Install MySQL

    - ``sudo apt install libmysqlclient-dev``
    - ``sudo apt install mysql-server``


#. Install R 4.0 . If you have other instances of R, we recommend uninstalling them.

    - ``sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9``
    - Add R's repository:

      - For 18.04 do: ``sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'``
      - For 20.04 do: ``sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'``
    - ``sudo apt update``
    - ``sudo apt install r-base``

#. Install Pandoc

    - ``sudo apt install pandoc``

#. Install GIT

    - ``sudo apt install git``

#. Install miniconda using these instructions_ 

#. Restart your current shell

#. Clone our repo:

    - ``git clone https://github.com/carissalow/rapids``

#. Create a python virtual environment:

    - ``cd rapids``
    - ``conda env create -f environment.yml -n MY_ENV_NAME``
    - ``conda activate MY_ENV_NAME``

#. Install R packages and virtual environment:

    - ``snakemake -j1 renv_install``
    - ``snakemake -j1 renv_restore``

      - This step could take several minutes to complete, especially if you have less than 3Gb of RAM or packages need to be compiled from source. Please be patient and let it run until completion. 

#. See Usage section below.


.. _usage-section:

Usage
======
Once RAPIDS is installed, follow these steps to start processing mobile data.

.. _db-configuration:

#. Configure the database connection:

   - Create an empty file called `.env` in the root directory (``rapids/``)
   - Add the following lines and replace your database-specific credentials (user, password, host, and database):

     .. code-block:: bash
        
         [MY_GROUP]
         user=MY_USER
         password=MY_PASSWORD
         host=MY_HOST
         port=3306
         database=MY_DATABASE

     .. note::

         ``MY_GROUP`` is a custom label for your credentials. It has to match ``DATABASE_GROUP`` in the ``config.yaml`` file_. It is not related to your database configuration.

#. Setup the participants' devices whose data you want to analyze, for this you have two options:

   #. **Automatically**. You can automatically include all devices that are stored in the ``aware_device`` table. If you want to control what devices and dates are included, see the Manual configuration::

        snakemake -j1 download_participants

   #. **Manually**. Create one file per participant in the ``rapids/data/external/`` directory. The file should NOT have an extension (i.e., no .txt). The name of the file will become the label for that participant in the pipeline.

      - The first line of the file should be the Aware ``device_id`` for that participant. If one participant has multiple device_ids (i.e. Aware had to be re-installed), add all device_ids separated by commas. 
      - The second line should list the device's operating system (``android`` or ``ios``). If a participant used more than one device (i.e., the participant changed phones and/or platforms mid-study) you can a) list each platform matching the order of the first line (``android,ios``), b) use ``android`` or ``ios`` if all phones belong to the same platform, or c) if you have an ``aware_device`` table in your database, set this line to ``multiple`` and RAPIDS will infer the multiple platforms automatically.
      - The third line is an optional human-friendly label that will appear in any plots for that participant.
      - The fourth line is optional and contains a start and end date separated by a comma ``YYYYMMDD,YYYYMMDD`` (e.g., ``20201301,20202505``). If these dates are specified, only data within this range will be processed, otherwise, all data from the device(s) will be used.

      For example, let's say participant `p01` had two AWARE device_ids and they were running Android between February 1st 2020 and March 3rd 2020. Their participant file would be named ``p01`` and contain:

        .. code-block:: bash

            3a7b0d0a-a9ce-4059-ab98-93a7b189da8a,44f20139-50cc-4b13-bdde-0d5a3889e8f9
            android
            Participant01
            2020/02/01,2020/03/03

#. Choose what features to extract:

   - See :ref:`Minimal Working Example<minimal-working-example>`.

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
