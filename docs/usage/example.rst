.. _analysis-workflow-example

Analysis Workflow Example
==========================

This is a quick guide for creating and running a simple pipeline to analysis an example dataset with 2 participants.

#. Install RAPIDS. See :ref:`Installation Section <install-page>`.

#. Make sure your database connection credentials in ``.env`` have write permission and set the correct ``MY_GROUP`` parameter in ``config.yaml`` file. See step 1 of :ref:`Usage Section <db-configuration>`.

#. Make sure your Conda (python) environment is active. See step 6 of :ref:`install-page`.

#. Run the following command to restore database from ``rapids_example.sql`` file. 

::

    sankemake -j1 restore_sql_file


#. Create example participants files with the following command.

::

    snakemake -j1 create_example_participant_files

#. Run the following command to analysis the example dataset.

    - Execute over a single core::

        snakemake -j1 --profile example_profile

    - Execute over multiple cores (here, we use 8 cores)::
    
        snakemake -j8 --profile example_profile

