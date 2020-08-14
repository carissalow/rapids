.. _analysis-workflow-example:

Analysis Workflow Example
==========================

This is a quick guide for creating and running a simple pipeline to analysis an example dataset with 2 participants.

#. Install RAPIDS. See :ref:`Installation Section <install-page>`.

#. Configure your database credentials (see the example below or step 1 of :ref:`Usage Section <db-configuration>` for more information).

    - Create an ``.env`` file at the root of RAPIDS folder 
    - Your MySQL user must have write permissions because we will restore our example database
    - Name your credentials group ``MY_GROUP``. 
    - If you are trying to connect to a local MySQL server from our docker container set your host according to this link_.
    - You can name your database any way you want, for example ``rapids_example``
    
    .. code-block::

        [MY_GROUP]
        user=rapids
        password=rapids
        host=127.0.0.1 # or use host.docker.internal from our docker container
        port=3306
        database=rapids_example

#. Make sure your conda environment is active (the environment is already active in our docker container). See step 6 of :ref:`install-page`.

#. Run the following command to restore database from ``rapids_example.sql`` file::

    snakemake -j1 restore_sql_file

#. Create example participants files with the following command::

    snakemake -j1 create_example_participant_files

#. Run the following command to analysis the example dataset.

    - Execute over a single core::

        snakemake -j1 --profile example_profile

    - Execute over multiple cores (here, we use 8 cores)::
    
        snakemake -j8 --profile example_profile

.. _link: https://stackoverflow.com/questions/24319662/from-inside-of-a-docker-container-how-do-i-connect-to-the-localhost-of-the-mach
