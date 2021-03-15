# Common Errors

##  Cannot connect to your MySQL server

???+ failure "Problem"
    ```bash
    **Error in .local(drv, \...) :** **Failed to connect to database: Error:
    Can\'t initialize character set unknown (path: compiled\_in)** :

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
    ```

???+ done "Solution"
    Please make sure the `DATABASE_GROUP` in `config.yaml` matches your DB credentials group in `.env`.

---

## Cannot start mysql in linux via `brew services start mysql`

???+ failure "Problem"
    Cannot start mysql in linux via `brew services start mysql`

???+ done "Solution"
    Use `mysql.server start`

---

## Every time I run force the download_dataset rule all rules are executed

???+ failure "Problem"
    When running `snakemake -j1 -R pull_phone_data` or `./rapids -j1 -R pull_phone_data` all the rules and files are re-computed

???+ done "Solution"
    This is expected behavior. The advantage of using `snakemake` under the hood is that every time a file containing data is modified every rule that depends on that file will be re-executed to update their results. In this case, since `download_dataset` updates all the raw data, and you are forcing the rule with the flag `-R` every single rule that depends on those raw files will be executed.
---

## Error `Table XXX doesn't exist` while running the `download_phone_data` or `download_fitbit_data` rule.

???+ failure "Problem"
    ```bash
    Error in .local(conn, statement, ...) : 
      could not run statement: Table 'db_name.table_name' doesn't exist
    Calls: colnames ... .local -> dbSendQuery -> dbSendQuery -> .local -> .Call
    Execution halted
    ```

???+ done "Solution"
    Please make sure the sensors listed in `[PHONE_VALID_SENSED_BINS][PHONE_SENSORS]` and the `[CONTAINER]` of each sensor you activated in `config.yaml`  match your database tables or files.

---
## How do I install RAPIDS on Ubuntu 16.04

???+ done "Solution"
    1.  Install dependencies (Homebrew - if not installed):
        -   `sudo apt-get install libmariadb-client-lgpl-dev libxml2-dev libssl-dev`
        -   Install [brew](https://docs.brew.sh/Homebrew-on-Linux) for linux and add the following line to `~/.bashrc`: `export PATH=$HOME/.linuxbrew/bin:$PATH`
        -   `source ~/.bashrc`

    1.  Install MySQL
        -   `brew install mysql`
        -   `brew services start mysql`

    2.  Install R, pandoc and rmarkdown:
        -   `brew install r`
        -   `brew install gcc@6` (needed due to this [bug](https://github.com/Homebrew/linuxbrew-core/issues/17812))
        -   `HOMEBREW_CC=gcc-6 brew install pandoc`

    3.  Install miniconda using these [instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

    4.  Clone our repo:
        -   `git clone https://github.com/carissalow/rapids`

    5.  Create a python virtual environment:
        -   `cd rapids`
        -   `conda env create -f environment.yml -n MY_ENV_NAME`
        -   `conda activate MY_ENV_NAME`

    6.  Install R packages and virtual environment:
        -   `snakemake renv_install`
        -   `snakemake renv_init`
        -   `snakemake renv_restore`

        This step could take several minutes to complete. Please be patient and let it run until completion.
---

## `mysql.h` cannot be found

???+ failure "Problem"
    ```bash
    --------------------------[ ERROR MESSAGE ]----------------------------
    <stdin>:1:10: fatal error: mysql.h: No such file or directory
    compilation terminated.
    -----------------------------------------------------------------------
    ERROR: configuration failed for package 'RMySQL'
    ```

???+ done "Solution"
    ```bash
    sudo apt install libmariadbclient-dev
    ```

---
## No package `libcurl` found

???+ failure "Problem"
    `libcurl` cannot be found

???+ done "Solution"
    Install `libcurl`
    ```bash
    sudo apt install libcurl4-openssl-dev
    ```

---
## Configuration failed because `openssl` was not found.

???+ failure "Problem"
    `openssl` cannot be found

???+ done "Solution"
    Install `openssl`
    ```bash
    sudo apt install libssl-dev
    ```
---
## Configuration failed because `libxml-2.0` was not found

???+ failure "Problem"
    `libxml-2.0` cannot be found

???+ done "Solution"
    Install `libxml-2.0`
    ```bash
    sudo apt install libxml2-dev
    ```

---
## SSL connection error when running RAPIDS

???+ failure "Problem"
    You are getting the following error message when running RAPIDS:
    ```bash
    Error: Failed to connect: SSL connection error: error:1425F102:SSL routines:ssl_choose_client_version:unsupported protocol.
    ```

???+ done "Solution"
    This is a bug in Ubuntu 20.04 when trying to connect to an old MySQL server with MySQL client 8.0. You should get the same error message if you try to connect from the command line. There you can add the option `--ssl-mode=DISABLED` but we can\'t do this from the R connector.

    If you can\'t update your server, the quickest solution would be to import your database to another server or to a local environment. Alternatively, you could replace `mysql-client` and `libmysqlclient-dev` with `mariadb-client` and `libmariadbclient-dev` and reinstall renv. More info about this issue [here](https://bugs.launchpad.net/ubuntu/+source/mysql-8.0/+bug/1872541)

---
## `DB_TABLES` key not found

???+ failure "Problem"
    If you get the following error `KeyError in line 43 of preprocessing.smk: 'PHONE_SENSORS'`, it means that the indentation of the key `[PHONE_SENSORS]` is not matching the other child elements of `PHONE_VALID_SENSED_BINS`
    
???+ done "Solution"
    You need to add or remove any leading whitespaces as needed on that line.

    ```yaml
    PHONE_VALID_SENSED_BINS:
        COMPUTE: False # This flag is automatically ignored (set to True) if you are extracting PHONE_VALID_SENSED_DAYS or screen or Barnett's location features
        BIN_SIZE: &bin_size 5 # (in minutes)
        PHONE_SENSORS: []
    ```

---
## Error while updating your conda environment in Ubuntu

???+ failure "Problem"
    You get the following error:
    ```bash
    CondaMultiError: CondaVerificationError: The package for tk located at /home/ubuntu/miniconda2/pkgs/tk-8.6.9-hed695b0_1003
        appears to be corrupted. The path 'include/mysqlStubs.h'
        specified in the package manifest cannot be found.
    ClobberError: This transaction has incompatible packages due to a shared path.
        packages: conda-forge/linux-64::llvm-openmp-10.0.0-hc9558a2_0, anaconda/linux-64::intel-openmp-2019.4-243
        path: 'lib/libiomp5.so'
    ```

???+ done "Solution"
    Reinstall conda

## Embedded nul in string

???+ failure "Problem"
    You get the following error when downloading sensor data:
    ```bash
    Error in result_fetch(res@ptr, n = n) : 
      embedded nul in string:
    ```

???+ done "Solution"
    This problem is due to the way `RMariaDB` handles a mismatch between data types in R and MySQL (see [this issue](https://github.com/r-dbi/RMariaDB/issues/121)). Since it seems this problem won't be handled by `RMariaDB`, you have two options:
    
    1. Remove the the null character from the conflictive table cell(s). You can adapt the following query on a MySQL server 8.0 or older
        ```sql
        update YOUR_TABLE set YOUR_COLUMN = regexp_replace(YOUR_COLUMN, '\0', '');
        ```
    2. If it's not feasible to modify your data you can try swapping `RMariaDB` with `RMySQL`. Just have in mind you might have problems connecting to modern MySQL servers running in Linux:
        - Add `RMySQL` to the renv environment by running the following command in a terminal open on RAPIDS root folder
        ```bash
        R -e 'renv::install("RMySQL")'
        ```
        - Go to `src/data/streams/pull_phone_data.R` or `src/data/streams/pull_fitbit_data.R` and replace `library(RMariaDB)` with `library(RMySQL)`
        - In the same file(s) replace `dbEngine <- dbConnect(MariaDB(), default.file = "./.env", group = group)` with `dbEngine <- dbConnect(MySQL(), default.file = "./.env", group = group)`
## There is no package called `RMariaDB`

???+ failure "Problem"
    You get the following error when executing RAPIDS:
    ```bash
    Error in library(RMariaDB) : there is no package called 'RMariaDB'
    Execution halted
    ```

???+ done "Solution"
    In RAPIDS v0.1.0 we replaced `RMySQL` R package with `RMariaDB`, this error means your R virtual environment is out of date, to update it run `snakemake -j1 renv_restore`
    
## Unrecognized output timezone "America/New_York"

???+ failure "Problem"
    When running RAPIDS with R 4.0.3 on MacOS on M1, lubridate may throw an error associated with the timezone.
    ```bash
    Error in C_force_tz(time, tz = tzone, roll):
       CCTZ: Unrecognized output timezone: "America/New_York"
    Calls: get_timestamp_filter ... .parse_date_time -> .strptime -> force_tz -> C_force_tz
    ```
???+ done "Solution"
   This is because R timezone library is not set. Please add `Sys.setenv(“TZDIR” = file.path(R.home(), “share”, “zoneinfo”))` to the file active.R in renv folder to set the timezone library. For further details on how to test if `TZDIR` is properly set, please refer to `https://github.com/tidyverse/lubridate/issues/928#issuecomment-720059233`. 
   
## Unimplemented MAX_NO_FIELD_TYPES

???+ failure "Problem"
    You get the following error when downloading Fitbit data:
    ```bash
    Error: Unimplemented MAX_NO_FIELD_TYPES
    Execution halted
    ```
???+ done "Solution"
    At the moment RMariaDB [cannot handle](https://github.com/r-dbi/RMariaDB/issues/127) MySQL columns of JSON type. Change the type of your Fitbit data column to `longtext` (note that the content will not change and will still be a JSON object just interpreted as a string).
    
## Running RAPIDS on Apple Silicon M1 Mac

???+ failure "Problem"
     You get the following error when installing pandoc or running rapids:
     ```bash
     MoSHI/rapids/renv/staging/1/00LOCK-KernSmooth/00new/KernSmooth/libs/KernSmooth.so: mach-0, but wrong architecture
     ```
???+ done "Solution"
    As of Feb 2020 in M1 macs, R needs to be installed via brew under Rosetta (x86 arch) due to some incompatibility with selected R libraries. To do this, run your terminal [via Rosetta](https://www.youtube.com/watch?v=nv2ylxro7rM&t=138s), then proceed with the usual brew installation command. x86 homebrew should be installed in `/usr/local/bin/brew `, you can check which brew you are using by typing `which brew`. Then use x86 homebrew to install R and restore RAPIDS packages (`renv_restore`). 
