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
