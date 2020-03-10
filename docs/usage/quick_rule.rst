Minimal Working Example 
=======================

The following is a quick guide for creating and running a simple pipeline to extract Call metrics for daily and night epochs of one participant monitored on the US East coast.

#. Make sure your database connection credentials in ``.env`` are correct. See step 1 of :ref:`Usage Section <db-configuration>`.

#. Create at least one participant file ``p01`` under ``data/external/``. See step 2 of :ref:`Usage Section <db-configuration>`.

#. Make sure your Conda (python) environment is active. See step 6 of :ref:`install-page`.

#. Replace the contents of the ``Snakefile`` with the following snippet
    
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


#. Modify the following settings in the ``config.yaml`` file with the values shown below (leave all other settings as they are)

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

#. Run the following command to execute RAPIDS

    ::

        snakemake

#. Daily and night call metrics will be found in files under the ``data/processed/p01/`` directory.

