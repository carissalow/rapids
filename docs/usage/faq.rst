Frequently Asked Questions 
============================

1. Cannot connect to the MySQL server
"""""""""""""""""""""""""""""""""""""""
**Error in .local(drv, ...) :**
**Failed to connect to database: Error: Can't initialize character set unknown (path: compiled_in)**
::

    Calls: dbConnect -> dbConnect -> .local -> .Call
    Execution halted
    [Tue Mar 10 19:40:15 2020]
    Error in rule download_dataset:
        jobid: 531
        output: data/raw/p60/locations_raw.csv

    RuleException:
    CalledProcessError in line 20 of /home/ubuntu/rapids/rules/preprocessing.snakefile:
    Command 'set -euo pipefail;  Rscript --vanilla /home/ubuntu/rapids/.snakemake/scripts/tmp_2jnvqs7.download_dataset.R' returned non-zero exit status 1.
    File "/home/ubuntu/rapids/rules/preprocessing.snakefile", line 20, in __rule_download_dataset
    File "/home/ubuntu/anaconda3/envs/moshi-env/lib/python3.7/concurrent/futures/thread.py", line 57, in run
    Shutting down, this might take some time.
    Exiting because a job execution failed. Look above for error message

**Solution:**

Please make sure the ``DATABASE_GROUP`` in ``config.yaml`` matches your DB credentials group in ``.env``.



2. Cannot start mysql in linux via ``brew services start mysql``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Use the following command instead:

``mysql.server start``


3. Every time I run ``snakemake -R download_dataset`` all rules are executed
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This is expected behavior. The advantage of using ``snakemake`` under the hood is that every time a file containing data is modified every rule that depends on that file will be re-executed to update their results. In this case, since ``download_dataset`` updates all the raw data, and you are forcing the rule with the flag ``-R`` every single rule that depends on those raw files will be executed.


4. Got an error ``Table XXX doesn't exist`` while running the download_dataset rule.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
::

    Error in .local(conn, statement, ...) : 
      could not run statement: Table 'db_name.table_name' doesn't exist
    Calls: colnames ... .local -> dbSendQuery -> dbSendQuery -> .local -> .Call
    Execution halted

**Solution:**
Please make sure the sensors listed in ``[PHONE_VALID_SENSED_BINS][TABLES]`` and each sensor section you activated in ``config.yaml`` match your database tables.



5. How do I install on Ubuntu 16.04
""""""""""""""""""""""""""""""""""""

#. Install dependencies (Homebrew - if not installed):

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
        - This step could take several minutes to complete. Please be patient and let it run until completion. 

#. See :ref:`Usage section <usage-section>`.



6. Configuration failed for package ``RMySQL``
""""""""""""""""""""""""""""""""""""""""""""""""
::

    --------------------------[ ERROR MESSAGE ]----------------------------
    <stdin>:1:10: fatal error: mysql.h: No such file or directory
    compilation terminated.
    -----------------------------------------------------------------------
    ERROR: configuration failed for package 'RMySQL'

Run ``sudo apt install libmariadbclient-dev``



7. No package ``libcurl`` found
"""""""""""""""""""""""""""""""""

The  ``libcurl`` needs to installed using the following command

Run ``sudo apt install libcurl4-openssl-dev``



8. Configuration failed because ``openssl`` was not found.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Install the ``openssl`` library using the following command

Run ``sudo apt install libssl-dev``


9. Configuration failed because ``libxml-2.0`` was not found
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Install the ``xml`` library using the following command 

Run ``sudo apt install libxml2-dev``

10. SSL connection error when running RAPIDS
""""""""""""""""""""""""""""""""""""""""""""""

You are getting the following error message when running RAPIDS:

``Error: Failed to connect: SSL connection error: error:1425F102:SSL routines:ssl_choose_client_version:unsupported protocol``.

This is a bug in Ubuntu 20.04 when trying to connect to an old MySQL server with MySQL client 8.0. You should get the same error message if you try to connect from the command line. There you can add the option ``--ssl-mode=DISABLED`` but we can't do this from the R connector.

If you can't update your server, the quickest solution would be to import your database to another server or to a local environment. Alternatively, you could replace ``mysql-client`` and ``libmysqlclient-dev`` with ``mariadb-client`` and ``libmariadbclient-dev`` and reinstall renv. More info about this issue here https://bugs.launchpad.net/ubuntu/+source/mysql-8.0/+bug/1872541

11. ``DB_TABLES`` key not found

If you get the following error ``KeyError in line 43 of preprocessing.smk: 'DB_TABLES'``, means that the indentation of the key ``DB_TABLES`` is not matching the other child elements of ``PHONE_VALID_SENSED_BINS`` and you need to add or remove any leading whitespaces as needed.

::

    PHONE_VALID_SENSED_BINS:
        COMPUTE: False # This flag is automatically ignored (set to True) if you are extracting PHONE_VALID_SENSED_DAYS or screen or Barnett's location features
        BIN_SIZE: &bin_size 5 # (in minutes)
        # Add as many sensor tables as you have, they all improve the computation of PHONE_VALID_SENSED_BINS and PHONE_VALID_SENSED_DAYS. 
        # If you are extracting screen or Barnett's location features, screen and locations tables are mandatory.
        DB_TABLES: []

12. Error while updating your conda environment in Ubuntu

If you get the following error try reinstalling conda

::
    CondaMultiError: CondaVerificationError: The package for tk located at /home/ubuntu/miniconda2/pkgs/tk-8.6.9-hed695b0_1003
        appears to be corrupted. The path 'include/mysqlStubs.h'
        specified in the package manifest cannot be found.
    ClobberError: This transaction has incompatible packages due to a shared path.
        packages: conda-forge/linux-64::llvm-openmp-10.0.0-hc9558a2_0, anaconda/linux-64::intel-openmp-2019.4-243
        path: 'lib/libiomp5.so'


.. ------------------------ Links --------------------------- ..

.. _bug: https://github.com/Homebrew/linuxbrew-core/issues/17812
.. _instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
.. _brew: https://docs.brew.sh/Homebrew-on-Linux
