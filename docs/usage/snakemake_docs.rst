.. _rapids-structure:

RAPIDS Structure
=================

.. _the-snakefile-file:

The ``Snakefile`` File
----------------------
The ``Snakefile`` file (see the actual `Snakefile`_) pulls the entire system together and can be thought of as the menu of RAPIDS allowing the user to define the sensor data that is desired. The first line in this file identifies the configuration file. Next are a list of included files that define the rules used to pull, clean, process, analyze and report on the data. Next is the ``all`` rule that list the sensor data (menu items) that would be processed by the pipeline. 

.. _includes-section:

Includes
"""""""""
There are 5 included files in the ``Snakefile`` file. 

    - ``renv.snakefile`` - This file defines the rules to manager the R packages that are used by RAPIDS. (See `renv`_)
    - ``preprocessing.snakefile`` - This file contains the rules that are used to preprocess the data such as downloading, cleaning and formatting. (See `preprocessing`_)
    - ``features.snakefile`` - This file contains the rules that used for behavioral feature extraction. (See `features`_)
    - ``models.snakefile`` - This file contains the rules that are used to build models from features that have been extreacted from the sensor data. (See `models`_)
    - ``reports.snakefile`` - The file contains the rules that are used to produce the reports based on the models produced. (See `reports`_)
    - ``mystudy.snakefile`` - The file contains the rules that you add that are specifically tailored to your project/study. (See `mystudy`_)

..  - ``analysis.snakefile`` - The rules that define how the data is analyzed is outlined in this file. (see `analysis <https://github.com/carissalow/rapids/blob/master/rules/analysis.snakefile>`_)
    
Includes are relative to the directory of the Snakefile in which they occur. For example, if above Snakefile resides in the directory ``my/dir``, then Snakemake will search for the include file at ``my/dir/path/to/other/snakefile``, regardless of the working directory.

.. _rule-all-section:

``Rule all:``
"""""""""""""
In RAPIDS the ``all`` rule indirectly specifies the features/sensors that are desired by listing the output files of the pipeline using the ``expand`` directive. The ``expand`` function allows the combination of different variables. Consider the following::

    expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),

If ``pids = ['p01','p02']`` and ``sensor = ['sms', 'calls']`` then the above directive would produce::

    ["data/raw/p01/sms_raw.csv", "data/raw/p01/calls_raw.csv", "data/raw/p02/sms_raw.csv", "data/raw/p02/calls_raw.csv"]

Thus, this allows the user of RAPIDS to define all of the desired output files without having to manually list all for the participants of the research. The way Snakemake works is that it looks for the rule that produces the desired output files and then executes that rule. For more information on ``expand`` see `The Expand Function`_


.. _the-env-file:

The ``.env`` File
-------------------
The database credentials for database server is placed in the .env file (Remember step 9 on :ref:`install-page` page). The format of the configurations are shown below::

    [MY_GROUP_NAME]
    user=MyUSER
    password=MyPassword
    host=MyIP
    port=3306


.. _the-config-file:

The ``config.yaml`` File
------------------------

The configurations for the pipeline are defined in the ``config.yaml`` (See `config.yaml`_). This contains global settings and variables that are used by the rules. Some of the global variables defined in the ``config.yaml`` file are briefly explained below:

    - ``SENSORS`` - This is a global variable that contains a list of the sensor/feature tables in the database that will be analyzed.
    - ``PIDS`` - This is the list of the participant IDs to include in the analysis. Create a file for each participant with a matching name ``pXXX`` containing the device_id in the ``data/external/`` directory. (Remember step 8 on the :ref:`install-page` page)
    - ``DAY_SEGMENTS`` - A variable used to list all of the common day segments. 
    - ``TIMEZONE`` - Time variable. Use timezone names from the `List of Timezone`_ and double check your code, for example EST is not US Eastern Time.
    - ``DATABASE_GROUP`` - Label for the database credentials group. (See :ref:`Configure the database connection <db-configuration>`.)
    - ``DOWNLOAD_DATASET`` - Variable used to store the name of the dataset that will be download for analysis. 

There are a number of other settings that are specific to the sensor/feature that will be pulled and analyzed by the pipeline. An example of the configuration settings for the :ref:`sms-sensor-doc` data is shown below::

    SMS:
        TYPES : [received, sent]
        FEATURES: 
            received: [count, distinctcontacts, timefirstsms, timelastsms, countmostfrequentcontact]
            sent: [count, distinctcontacts, timefirstsms, timelastsms, countmostfrequentcontact]
        DAY_SEGMENTS: *day_segments  

The ``TYPES`` setting defines the type of SMS data that will be analyzed. ``FEATURES`` defines the features of the data for each the type of SMS data being analyzed. Finally, ``DAY_SEGMENTS`` list the day segment (times of day) that the data is captured.

.. _rules-syntax:

The ``Rules`` Directory 
------------------------

The ``rules`` directory contains the ``snakefiles`` that were included in the ``Snakefile`` file. A short description of these files are given in the :ref:`includes-section` section. 


Rules
""""""

A Snakemake workflow is defined by specifying rules in a ``Snakefile`` (See the features_ snakefile as an actual example). Rules decompose the workflow into small steps (e.g., the application of a single tool) by specifying how to create sets of output files from sets of input files. Snakemake automatically determines the dependencies between the rules by matching file names. Thus, a rule can consist of a name, input files, output files, and a command to generate the output from the input. The following is the basic structure of a Snakemake rule::

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", "path/to/another/outputfile"
        script: "path/to/somescript.R"


A sample rule from the RAPIDS source code is shown below::

    rule sms_features:
        input: 
            "data/raw/{pid}/messages_with_datetime.csv"
        params:
            sms_type = "{sms_type}",
            day_segment = "{day_segment}",
            features = lambda wildcards: config["SMS"]["FEATURES"][wildcards.sms_type]
        output:
            "data/processed/{pid}/sms_{sms_type}_{day_segment}.csv"
        script:
            "../src/features/sms_features.R"


The ``rule`` directive specifies the name of the rule that is being defined. ``params`` defines the additional parameters that needs to be set for the rule. In the example immediately above, the parameters will be pasted to the script defined in the ``script`` directive of the rule. Instead of ``script`` a ``shell`` command call can also be called by replacing the ``script`` directive of the rule and replacing it with the lines similar to the folllowing::

        shell: "somecommand {input} {output}"

Here input and output (and in general any list or tuple) automatically evaluate to a space-separated list of files (i.e. ``path/to/inputfile path/to/other/inputfile``).  It should be noted that rules can defined without input and output as seen in the ``renv.snakemake``. For more information see `Rules documentation`_ and for an actual example see the `renv`_ snakefile.

.. _wildcards:

Wildcards
""""""""""
There are times that it would be useful to generalize a rule to be applicable to a number of e.g. datasets. For this purpose, wildcards can be used. Consider the sample code from above again repeated below for quick reference.::

    rule sms_features:
        input: 
            "data/raw/{pid}/messages_with_datetime.csv"
        params:
            sms_type = "{sms_type}",
            day_segment = "{day_segment}",
            features = lambda wildcards: config["SMS"]["FEATURES"][wildcards.sms_type]
        output:
            "data/processed/{pid}/sms_{sms_type}_{day_segment}.csv"
        script:
            "../src/features/sms_features.R"

If the rule’s output matches a requested file, the substrings matched by the wildcards are propagated to the input and params directives. For example, if another rule in the workflow requires the file ``data/processed/p01/sms_sent_daily.csv``, Snakemake recognizes that the above rule is able to produce it by setting ``pid=p01``, ``sms_type=sent`` and ``day_segment=daily``. Thus, it requests the input file ``data/raw/p01/messages_with_datetime.csv`` as input, sets ``sms_type=sent``, ``day_segment=daily`` in the ``params`` directive and executes the script. ``../src/features/sms_features.R``. See the preprocessing_ snakefile for an actual example. 


.. _the-data-directory:

The ``data`` Directory
-----------------------

This directory contains the data files for the project. These directories are as follows:

    - ``external`` - This directory stores the participant `pxxx` files that contains the device_id and the type of device as well as data from third party sources. (Remember step 8 on :ref:`install-page` page)
    - ``raw`` - This directory contains the original, immutable data dump from the sensor database.
    - ``interim`` - This directory would contain intermediate data that has been transformed but has not been completely analyzed.
    - ``processed`` - This directory contains the final canonical data sets for modeling.


.. _the-src-directory:

The ``src`` Directory
----------------------

The ``src`` directory holds all of the scripts used by the pipeline for data manipulation. These scripts can be in any programming language including but not limited to Python_, R_ and Julia_. This directory is organized into the following directories:

    - ``data`` - This directory contains scripts that are used to download and preprocess raw data that will be used in analysis. See `data directory`_
    - ``features`` - This directory contains scripts to extract behavioral features. See `features directory`_
    - ``models`` - This directory contains the model scripts for building and training models. See `models directory`_
    - ``visualization`` - This directory contains the scripts to create plots and reports that visualize the results of the models. See `visualization directory`_


.. _the-report-directory:

The ``reports`` Directory
--------------------------

This contains the reports of the results of the analysis done by the pipeline. 

    .. _Python: https://www.python.org/
    .. _Julia: https://julialang.org/
    .. _R: https://www.r-project.org/
    .. _`List of Timezone`: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
    .. _`The Expand Function`: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-expand-function
    .. _`example snakefile`: https://github.com/carissalow/rapids/blob/master/rules/features.snakefile
    .. _renv: https://github.com/carissalow/rapids/blob/master/rules/renv.snakefile
    .. _preprocessing: https://github.com/carissalow/rapids/blob/master/rules/preprocessing.snakefile
    .. _features: https://github.com/carissalow/rapids/blob/master/rules/features.snakefile
    .. _models: https://github.com/carissalow/rapids/blob/master/rules/models.snakefile
    .. _reports: https://github.com/carissalow/rapids/blob/master/rules/reports.snakefile
    .. _mystudy: https://github.com/carissalow/rapids/blob/master/rules/mystudy.snakefile
    .. _`Rules documentation`: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rules
    .. _`data directory`: https://github.com/carissalow/rapids/tree/master/src/data
    .. _`features directory`: https://github.com/carissalow/rapids/tree/master/src/features
    .. _`models directory`: https://github.com/carissalow/rapids/tree/master/src/models
    .. _`visualization directory`: https://github.com/carissalow/rapids/tree/master/src/visualization
    .. _`config.yaml`: https://github.com/carissalow/rapids/blob/master/config.yaml
    .. _`Snakefile`: https://github.com/carissalow/rapids/blob/master/Snakefile


::

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── config.yaml        <- The configuration settings for the pipeline.
    ├── environment.yml    <- Environmental settings - channels and dependences that are installed in the env)
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── packrat            <- Installed R dependences. (Packrat is a dependency management system for R) 
    │                         (Depreciated - replaced by renv)
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── renv.lock          <- List of R packages and dependences for that are installed for the pipeline.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting.
    │
    ├── rules              
    │   ├── features       <- Rules to process the feature data pulled in to pipeline.
    │   ├── models         <- Rules for building models.
    │   ├── mystudy        <- Rules added by you that are specifically tailored to your project/study.
    │   ├── packrat        <- Rules for setting up packrat. (Depreciated replaced by renv)
    │   ├── preprocessing  <- Preprocessing rules to clean data before processing.
    │   ├── renv           <- Rules for setting up renv and R packages.
    │   └── reports        <- Snakefile used to produce reports.
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── Snakemake          <- The root snakemake file (the equivalent of a Makefile)
    ├── src                <- Source code for use in this project. Can be in any language e.g. Python, 
    │   │                     R, Julia, etc.
    │   │
    │   ├── data           <- Scripts to download or generate data. Can be in any language e.g. Python, 
    │   │                     R, Julia, etc.
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling. Can be in any language 
    │   │                     e.g. Python, R, Julia, etc.
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make prediction. Can 
    │   │                     be in any language e.g. Python, R, Julia, etc.
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations. Can be 
    │                         in any language e.g. Python, R, Julia, etc.
    ├── tests
    │   ├── data           <- Replication of the project root data directory for testing.
    │   ├── scripts        <- Scripts for testing.
    │   ├── settings       <- The config and settings files for running tests.
    │   └── Snakefile      <- The Snakefile for testing only.
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.testrun.org