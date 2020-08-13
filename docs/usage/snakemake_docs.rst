.. _rapids-structure:

RAPIDS Structure
=================

.. _the-config-file:

The ``config.yaml`` File
------------------------

RAPIDS configuration settings are defined in ``config.yaml`` (See `config.yaml`_). This is the only file that you need to understand in order to compute the features that RAPIDS ships with. 

It has global settings like ``PIDS``, ``DAY_SEGMENTS``, among others (see :ref:`global-sensor-doc` for more information). As well as per sensor settings, for example, for the :ref:`messages-sensor-doc`::

      | ``MESSAGES:``
      |     ``COMPUTE: True``
      |      ``DB_TABLE: messages``
      |      ``...``

.. _the-snakefile-file:

The ``Snakefile`` File
----------------------
The ``Snakefile`` file (see the actual `Snakefile`_) pulls the entire system together. The first line in this file identifies the configuration file. Next are a list of included directives that import the rules used to pull, clean, process, analyze and report data. It compiles the list of ``files_to_compute`` by scaning the config file looking for the sensors with a ``COMPUTE`` flag equal to ``True``. Then, the ``all`` rule is called with this list which prompts Snakemake to exectue the pipeline (raw files, intermediate files, feature files, reports, etc). 

.. _includes-section:

Includes
"""""""""
There are 5 included files in the ``Snakefile`` file. 

    - ``renv.smk`` - Rules to create, backup and restore the R renv virtual environment for RAPIDS. (See `renv`_)
    - ``preprocessing.smk`` - Rules that are used to pre-preprocess the data such as downloading, cleaning and formatting. (See `preprocessing`_)
    - ``features.smk`` - Rules that used for behavioral feature extraction. (See `features`_)
    - ``models.smk`` - Rules that are used to build models from features that have been extreacted from the sensor data. (See `models`_)
    - ``reports.smk`` - Rules that are used to produce reports and visualizations. (See `reports`_)
    
Includes are relative to the root directory.

.. _rule-all-section:

``Rule all:``
"""""""""""""
In RAPIDS the ``all`` rule lists the output files we expect the pipeline to compute. Before the ``all`` rule is called snakemake checks the ``config.yaml`` and adds all the rules for which the sensors ``COMPUTE`` parameter is ``True``. The ``expand`` function allows us to generate a list of file paths that have a common structure except for PIDS or other parameters. Consider the following::

    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["MESSAGES"]["DB_TABLE"]))

If ``pids = ['p01','p02']`` and ``config["MESSAGES"]["DB_TABLE"] = messages`` then the above directive would produce::

    ["data/raw/p01/messages_raw.csv", "data/raw/p02/messages_raw.csv"]

Thus, this allows us to define all the desired output files without having to manually list each path for every participant and every sensor. The way Snakemake works is that it looks for the rule that produces the desired output files and then executes that rule. For more information on ``expand`` see `The Expand Function`_


.. _the-env-file:

The ``.env`` File
-------------------
Your database credentials are stored in the ``.env`` file (See :ref:`install-page`)::

    [MY_GROUP_NAME]
    user=MyUSER
    password=MyPassword
    host=MyIP/DOMAIN
    port=3306

.. _rules-syntax:

The ``Rules`` Directory 
------------------------

The ``rules`` directory contains the ``snakefiles`` that were included in the main ``Snakefile`` file. A short description of these files are given in the :ref:`includes-section` section. 


Rules
""""""

A Snakemake workflow is defined by rules (See the features_ snakefile as an actual example). Rules decompose the workflow into small steps by specifying what output files should be created by running a script on a set of input files. Snakemake automatically determines the dependencies between the rules by matching file names. Thus, a rule can consist of a name, input files, output files, and a command to generate the output from the input. The following is the basic structure of a Snakemake rule::

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", "path/to/another/outputfile"
        script: "path/to/somescript.R"


A sample rule from the RAPIDS source code is shown below::

    rule messages_features:
        input: 
            expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["MESSAGES"]["DB_TABLE"])
        params:
            messages_type = "{messages_type}",
            day_segment = "{day_segment}",
            features = lambda wildcards: config["MESSAGES"]["FEATURES"][wildcards.messages_type]
        output:
            "data/processed/{pid}/messages_{messages_type}_{day_segment}.csv"
        script:
            "../src/features/messages_features.R"


The ``rule`` directive specifies the name of the rule that is being defined. ``params`` defines additional parameters for the rule's script. In the example above, the parameters are passed to the ``messages_features.R`` script as an dictionary. Instead of ``script`` a ``shell`` command call can also be called by replacing the ``script`` directive of the rule and replacing it with::

        shell: "somecommand {input} {output}"

It should be noted that rules can be defined without input and output as seen in the ``renv.snakemake``. For more information see `Rules documentation`_ and for an actual example see the `renv`_ snakefile.

.. _wildcards:

Wildcards
""""""""""
There are times when the same rule should be applied to different participants and day segments. For this we use wildcards ``{my_wildcard}``. All wildcards are inferred from the files listed in the ``all` rule of the ``Snakefile`` file and therefore from the output of any rule::

    rule messages_features:
        input: 
            expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["MESSAGES"]["DB_TABLE"])
        params:
            messages_type = "{messages_type}",
            day_segment = "{day_segment}",
            features = lambda wildcards: config["MESSAGES"]["FEATURES"][wildcards.messages_type]
        output:
            "data/processed/{pid}/messages_{messages_type}_{day_segment}.csv"
        script:
            "../src/features/messages_features.R"

If the rule’s output matches a requested file, the substrings matched by the wildcards are propagated to the input and params directives. For example, if another rule in the workflow requires the file ``data/processed/p01/messages_sent_daily.csv``, Snakemake recognizes that the above rule is able to produce it by setting ``pid=p01``, ``messages_type=sent`` and ``day_segment=daily``. Thus, it requests the input file ``data/raw/p01/messages_with_datetime.csv`` as input, sets ``messages_type=sent``, ``day_segment=daily`` in the ``params`` directive and executes the script. ``../src/features/messages_features.R``. See the preprocessing_ snakefile for an actual example. 


.. _the-data-directory:

The ``data`` Directory
-----------------------

This directory contains the data files for the project. These directories are as follows:

    - ``external`` - This directory stores the participant `pxxx` files as well as data from third party sources (see :ref:`install-page` page).
    - ``raw`` - This directory contains the original, immutable data dump from your database.
    - ``interim`` - This directory contains intermediate data that has been transformed but do not represent features.
    - ``processed`` - This directory contains all behavioral features.


.. _the-src-directory:

The ``src`` Directory
----------------------

The ``src`` directory holds all the scripts used by the pipeline for data manipulation. These scripts can be in any programming language including but not limited to Python_, R_ and Julia_. This directory is organized into the following directories:

    - ``data`` - This directory contains scripts that are used to download and preprocess raw data that will be used in analysis. See `data directory`_
    - ``features`` - This directory contains scripts to extract behavioral features. See `features directory`_
    - ``models`` - This directory contains the scripts for building and training models. See `models directory`_
    - ``visualization`` - This directory contains the scripts to create plots and reports. See `visualization directory`_


.. _RAPIDS_directory_structure:

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
