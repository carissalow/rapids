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

Please make sure the ``MY_GROUP`` in ``config.yaml`` and ``.env`` match.



2. Cannot start mysql in linux via ``brew services start mysql``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
use the following command instead:

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
Please make sure the sensors listed in SENSORS in ``config.yaml`` match your database tables.



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

sudo apt install libmariadbclient-dev



7. No package ``libcurl`` found
"""""""""""""""""""""""""""""""""

The  ``libcurl`` needs to installed using the following command

``sudo apt install libcurl4-openssl-dev``



8. Configuration failed because ``openssl`` was not found.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Install the ``openssl`` library using the following command

``sudo apt install libssl-dev``


9. Configuration failed because ``libxml-2.0`` was not found
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Install the ``xml`` library using the following command 

``sudo apt install libxml2-dev``

.. ------------------------ Links --------------------------- ..

.. _bug: https://github.com/Homebrew/linuxbrew-core/issues/17812
.. _instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
.. _brew: https://docs.brew.sh/Homebrew-on-Linux
