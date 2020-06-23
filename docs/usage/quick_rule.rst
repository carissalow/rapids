.. _minimal-working-example:

Minimal Working Example 
=======================

This is a quick guide for creating and running a simple pipeline to extract call features for daily and night epochs of one participant monitored on the US East coast.

#. Make sure your database connection credentials in ``.env`` are correct. See step 1 of :ref:`Usage Section <db-configuration>`.

#. Create at least one participant file ``p01`` under ``data/external/``. See step 2 of :ref:`Usage Section <db-configuration>`.

#. Make sure your Conda (python) environment is active. See step 6 of :ref:`install-page`.

#. Modify the following settings in the ``config.yaml`` file with the values shown below (leave all other settings as they are)

    ::
        PIDS: [p01]
        
        DAY_SEGMENTS: &day_segments
            [daily, night]

        TIMEZONE: &timezone
            America/New_York
        
        DATABASE_GROUP: &database_group
            MY_GROUP (change this if you added your DB credentials to .env with a different label)

        CALLS:
            COMPUTE: True
            DB_TABLE: calls (only change DB_TABLE if your database calls table has a different name)
    
    For more information on the ``calls`` sensor see :ref:`call-sensor-doc`

#. Run the following command to execute RAPIDS

    ::

        snakemake -j1

#. Daily and night call metrics will be found in files under the ``data/processed/p01/`` directory.

