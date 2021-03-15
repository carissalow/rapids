# Migration guides

## Migrating from RAPIDS 0.4.x or older

There are four actions that you need to take if you were using RAPIDS `0.4.3` or older ([before Feb 9th, 2021](https://github.com/carissalow/rapids/releases/tag/v0.4.3)):

??? check "Check the new Overview page"
      Check the new [Overview](../setup/overview/) page. Hopefully, it is a better overview of RAPIDS and provides answers to Frequently Asked Questions.

??? check "Deploy RAPIDS in a new folder"

      - Clone RAPIDS 1.x in a new folder (do not pull the updates in your current folder)
      - Activate your conda environment
      - Install renv again `snakemake -j1 renv_install` (for Ubuntu take advantage of the [platform specific R `renv` instructions](../setup/installation))
      - Restore renv packages `snakemake -j1 renv_restore` (for Ubuntu take advantage of the [platform specific R `renv` instructions](../setup/installation))
      - Move your participant files `pxx.yaml` to the new folder
      - Move your time segment files to the new folder
      - Move your `.env` file to the new folder

??? check "Migrate your `.env` file to the new `credentials.yaml` format"
      The `.env` file is not used anymore, the same credential groups are stored in `credentials.yaml`, migrate your `.env` file by running:
      ```bash
      python tools/update_format_env.py
      ```

??? check "Reconfigure your `config.yaml`"
      Reconfigure your `config.yaml` file by hand (don't copy and paste the old one). Some keys and values changed but the defaults should be compatible with the things you know from RAPIDS 0.x (see below).

The most relevant changes to RAPIDS that you need to know about are:

??? danger "We introduced the concept of data streams"

      RAPIDS abstracts sensor data logged by different devices, platforms and stored in different data containers as [data streams](../datastreams/data-streams-introduction/).

      The default data stream for `PHONE` is [`aware_mysql`](../datastreams/aware-mysql/), and the default for `FITBIT` is [`fitbitjson_mysql`](../datastreams/fitbitjson-mysql/). This is compatible with the old functionality (AWARE and JSON Fitbit data stored in MySQL). These values are set in `[PHONE_DATA_STREAMS][USE]` and `[FITBIT_DATA_STREAMS][USE]`.

      You can [add new data stream](../datastreams/add-new-data-streams/) formats (sensing apps) and containers (database engines, file types, etc.).
      
      If you were processing your Fitbit data either in JSON or plain text (parsed) format, and it was stored in MySQL or CSV files, the changes that you made to your raw data will be compatible. Just choose [`fitbitjson_mysql`](../datastreams/fitbitjson-mysql/), [`fitbitparsed_mysql`](../datastreams/fitbitparsed-mysql/), [`fitbitjson_csv`](../datastreams/fitbitjson-csv/), [`fitbitparsed_csv`](../datastreams/fitbitparsed-csv/) accordingly and set it in `[FITBIT_DATA_STREAMS][USE]`. 
      
      In the future, you will not have to change your raw data; you will be able to just change column mappings/values in the data stream's `format.yaml` file.

??? danger "We introduced multiple time zones"
      You can now process data from participants that visited multiple time zones. The default is still a single time zone (America/New_York). See how to handle [multiple time zones](../setup/configuration/#multiple-timezones)

??? danger "The keyword `multiple` is now `infer`"
      When processing data from smartphones, RAPIDS allows you to [infer](../setup/configuration/#participant-files) the OS of a smartphone by using the keyword `multiple` in the `[PLATFORM]` key of participant files. Now RAPIDS uses `infer` instead of `multiple` Nonetheless, `multiple` still works for backward compatibility.

??? danger "A global `DATABASE_GROUP` does not exist anymore"
      There is no global `DATABASE_GROUP` anymore. Each data stream that needs credentials to connect to a database has its own [`DATABASE_GROUP` config key](../setup/configuration/#data-stream-configuration). The groups are defined in `credentials.yaml` instead of the `.env`.

??? danger "`[DEVICE_SENSOR][TABLE]` is now `[DEVICE_SENSOR][CONTAINER]`"
      We renamed the keys `[DEVICE_SENSOR][TABLE]` to `[DEVICE_SENSOR][CONTAINER]` to reflect that, with the introduction of data streams, they can point to a database table, file, or any other data container.

??? danger "Creating participant files from the AWARE_DEVICE_TABLE is deprecated"
    In previous versions of RAPIDS, you could create participant files automatically using the `aware_device` table. We deprecated this option but you can still achieve the same results if you export the output of the following SQL query as a CSV file and follow the instructions to [create participant files from CSV files](../setup/configuration/#automatic-creation-of-participant-files):
    
    ```sql
    SELECT device_id, device_id as fitbit_id, CONCAT("p", _id) as empatica_id, CONCAT("p", _id) as pid, if(brand = "iPhone", "ios", "android") as platform, CONCAT("p", _id)  as label, DATE_FORMAT(FROM_UNIXTIME((timestamp/1000)- 86400), "%Y-%m-%d") as start_date, CURRENT_DATE as end_date from aware_device order by _id;
    ```
??? danger "`SCR_SCRIPT` and `SRC_LANGUAGE` are replaced by `SRC_SCRIPT`"
    The attributes `SCR_SCRIPT` and `SRC_LANGUAGE` of every sensor `PROVIDER` are replaced by `SRC_SCRIPT`. `SRC_SCRIPT` is a relative path from the RAPIDS root folder to that provider's feature script. We did this to simplify and clarify where the features scripts are stored. 
    
    There are no actions to take unless you created your own feature provider; update it with your feature script path.
## Migrating from RAPIDS beta

If you were relying on the [old docs](https://rapidspitt.readthedocs.io/en/latest/) and the most recent version of RAPIDS you are working with is from or before [Oct 13, 2020](https://github.com/carissalow/rapids/commit/640890c7b49492d150accff5c87b1eb25bd97a49) you are using the beta version of RAPIDS.

You can start using the RAPIDS `0.1.0` right away, just take into account the following:

??? check "Deploy RAPIDS in a new folder"
      - [Install](setup/installation.md) a new copy of RAPIDS (the R and Python virtual environments didn't change so the cached versions will be reused)
      - Make sure you don't skip a new Installation step to give execution permissions to the RAPIDS script: `chmod +x rapids`
      - Move your old `.env` file
      - Move your participant files

??? check "Migrate your participant files"
      You can migrate your old participant files to the new YAML format:
      ```bash
      python tools/update_format_participant_files.py
      ```

??? check "Follow the new Configuration guide"
      Follow the new [Configuration](https://www.rapids.science/0.1/setup/configuration/) guide

??? check "Learn more about the new way to run RAPIDS"
      Get familiar with the new way of [Executing](https://www.rapids.science/0.1/setup/execution) RAPIDS
