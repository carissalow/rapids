Quick Rule 
=============

The following is a quick guide for creating and running a simple rule.

#. Setup database connection credential the ``.env``. See step 1 under the :ref:`Usage Section of Install <db-configuration>` page.

#. Create at least one participant file ``p01`` in the ``data/external``. See step 2 under the :ref:`Usage Section of Install <db-configuration>` page.

#. Activate the Conda (python) environment. See install step 2 on :ref:`install-page` page.

#. Prepare the ``Snakefile`` by replacing the code in the ``Snakefile`` with the following
    
    ::

        configfile: "config.yaml"
        include: "rules/packrat.snakefile"
        include: "rules/preprocessing.snakefile"
        include: "rules/features.snakefile"
        include: "rules/reports.snakefile"

        rule all:
            input:
                expand("data/processed/{pid}/call_{call_type}_{day_segment}.csv",
                                pid=config["PIDS"], 
                                call_type=config["CALLS"]["TYPES"],
                                day_segment = config["CALLS"]["DAY_SEGMENTS"]),


#. Prepare the ``config.yaml`` changing the following settings only to values shown below (leave all other settings as they are)

    ::

        SENSORS: [calls]

        FITBIT_TABLE: []
        FITBIT_SENSORS: []

        PIDS: [p01]
        
        DAY_SEGMENTS: &day_segments
            [daily, night]

        TIMEZONE: &timezone
            America/New_York
        
        DATABASE_GROUP: &database_group
            MY_GROUP
    
    For more information on the ``calls`` sensor see :ref:`call-sensor-doc`

#. Run the following command to execute the rule

    ::

        snakemake

#. The results of the execution of the rule will be found in the ``data/processed/p01`` directory.

