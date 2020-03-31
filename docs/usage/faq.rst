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

Please make sure the ``MY_GROUP`` in ``config.yaml`` and ``.env`` are the same.

2. Cannot start mysql for linux via ``brew services start mysql``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
use the following code instead:

``mysql.server start``

3. Every time I run ``snakemake -R download_dataset`` all rules are executed
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This is expected behavior. The advantage of using ``snakemake`` under the hood is that every time a file containing data is modified every rule that depends on that file will be re-executed to update their results. In this case, since ``download_dataset`` updates all the raw data, every single rule that depends on those raw files will be executed.

4. Got an error while running ``snakemake packrat_install`` to setup the RAPIDS environment
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
**Error:**
::

    SyntaxError in line 19 of /Users/caomz/rapids/Snakefile:
    Unexpected keyword expand in rule definition (Snakefile, line 19)

**Solution:**

Please make sure there are no extra whitespaces in Snakefile.
