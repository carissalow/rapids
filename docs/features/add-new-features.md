# Add New Features

!!! hint
    - We recommend reading the [Behavioral Features Introduction](../feature-introduction/) before reading this page.
    - You can implement new features in Python or R scripts.
    - You won't have to deal with time zones, dates, times, data cleaning, or preprocessing. The data that RAPIDS pipes to your feature extraction code are ready to process.

## New Features for Existing Sensors

You can add new features to any existing sensors (see list below) by adding a new provider in three steps:

1. [Modify](#modify-the-configyaml-file) the `config.yaml` file 
2. [Create](#create-a-feature-provider-script) your feature provider script
3. [Implement](#implement-your-feature-extraction-code) your features extraction code
   
As a tutorial, we will add a new provider for `PHONE_ACCELEROMETER` called `VEGA` that extracts `feature1`, `feature2`, `feature3` with a Python script that requires a parameter from the user called `MY_PARAMETER`.

??? info "Existing Sensors"
    An existing sensor of any device with a configuration entry in `config.yaml`:

    Smartphone (AWARE)

    - Phone Accelerometer
    - Phone Activity Recognition
    - Phone Applications Crashes
    - Phone Applications Foreground
    - Phone Applications Notifications
    - Phone Battery
    - Phone Bluetooth
    - Phone Calls
    - Phone Conversation
    - Phone Data Yield
    - Phone Keyboard
    - Phone Light
    - Phone Locations
    - Phone Log
    - Phone Messages
    - Phone Screen
    - Phone WiFI Connected
    - Phone WiFI Visible

    Fitbit

    - Fitbit Data Yield
    - Fitbit Heart Rate Summary
    - Fitbit Heart Rate Intraday
    - Fitbit Sleep Summary
    - Fitbit Sleep Intraday
    - Fitbit Steps Summary
    - Fitbit Steps Intraday

    Empatica

    - Empatica Accelerometer
    - Empatica Heart Rate
    - Empatica Temperature
    - Empatica Electrodermal Activity
    - Empatica Blood Volume Pulse
    - Empatica Inter Beat Interval
    - Empatica Tags


### Modify the `config.yaml` file

In this step, you need to add your provider configuration section under the relevant sensor in `config.yaml`. See our example for our tutorial's `VEGA` provider for  `PHONE_ACCELEROMETER`:

??? example "Example configuration for a new accelerometer provider `VEGA`"
    ```yaml hl_lines="12 13 14 15 16"
    PHONE_ACCELEROMETER:
        CONTAINER: accelerometer
        PROVIDERS:
            RAPIDS: # this is a feature provider
                COMPUTE: False
                ...
            
            PANDA: # this is another feature provider
                COMPUTE: False
                ...

            VEGA: # this is our new feature provider
                COMPUTE: False
                FEATURES: ["feature1", "feature2", "feature3"]
                MY_PARAMTER: a_string
                SRC_SCRIPT: src/features/phone_accelerometer/vega/main.py
            
    ```

| Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description
|---|---|
|`[COMPUTE]`| Flag to activate/deactivate your provider
|`[FEATURES]`| List of features your provider supports. Your provider code should only return the features on this list
|`[MY_PARAMTER]`| An arbitrary parameter that our example provider `VEGA` needs. This can be a boolean, integer, float, string, or an array of any of such types.
|`[SRC_SCRIPT]`| The relative path from RAPIDS' root folder to a script that computes the features for this provider. It can be implemented in R or Python.

### Create a feature provider script

Create your feature Python or R script called `main.py` or `main.R` in the correct folder, `src/feature/[sensorname]/[providername]/`. RAPIDS automatically loads and executes it based on the config key `[SRC_SCRIPT]` you added in the last step. For our example, this script is:
```bash
src/feature/phone_accelerometer/vega/main.py
```

### Implement your feature extraction code
Every feature script (`main.[py|R]`) needs a `[providername]_features` function with specific parameters. RAPIDS calls this function with the sensor data ready to process and with other functions and arguments you will need.

=== "Python function"
    ```python
    def [providername]_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):
        # empty for now
        return(your_features_df)
    ```

=== "R function"
    ```r
    [providername]_features <- function(sensor_data, time_segment, provider){
        # empty for now
        return(your_features_df)
    }
    ```

| Parameter&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description
|---|---|
|`sensor_data_files`| Path to the CSV file containing the data of a single participant. This data has been cleaned and preprocessed. Your function will be automatically called for each participant in your study (in the `[PIDS]` array in `config.yaml`) 
|`time_segment`| The label of the time segment that should be processed.
|`provider`| The parameters you configured for your provider in `config.yaml` will be available in this variable as a dictionary in Python or a list in R. In our example, this dictionary contains `{MY_PARAMETER:"a_string"}`
|`filter_data_by_segment`| Python only. A function that you will use to filter your data. In R, this function is already available in the environment.
|`*args`| Python only. Not used for now
|`**kwargs`| Python only. Not used for now


The next step is to implement the code that computes your behavioral features in your provider script's function. As with any other script, this function can call other auxiliary methods, but in general terms, it should have three stages:

??? info "1. Read a participant's data by loading the CSV data stored in the file pointed by `sensor_data_files`"
    ``` python
    acc_data = pd.read_csv(sensor_data_files["sensor_data"])
    ```

    Note that the phone's battery, screen, and activity recognition data are given as episodes instead of event rows (for example, start and end timestamps of the periods the phone screen was on)


??? info "2. Filter your data to process only those rows that belong to `time_segment`"

    This step is only one line of code, but keep reading to understand why we need it.
    ```python
    acc_data = filter_data_by_segment(acc_data, time_segment)
    ```

    You should use the `filter_data_by_segment()` function to process and group those rows that belong to each of the [time segments RAPIDS could be configured with](../../setup/configuration/#time-segments).

    Let's understand the `filter_data_by_segment()` function with an example. A RAPIDS user can extract features on any arbitrary [time segment](../../setup/configuration/#time-segments). A time segment is a period that has a label and one or more instances. For example, the user (or you) could have requested features on a daily, weekly, and weekend basis for `p01`. The labels are arbitrary, and the instances depend on the days a participant was monitored for: 

     - the daily segment could be named `my_days` and if `p01` was monitored for 14 days, it would have 14 instances
     - the weekly segment could be named `my_weeks` and if `p01` was monitored for 14 days, it would have 2 instances.
     - the weekend segment could be named `my_weekends` and if `p01` was monitored for 14 days, it would have 2 instances.
    
    For this example, RAPIDS will call your provider function three times for `p01`, once where `time_segment` is `my_days`, once where `time_segment` is `my_weeks`, and once where `time_segment` is `my_weekends`. In this example, not every row in `p01`'s data needs to take part in the feature computation for either segment **and** the rows need to be grouped differently. 
    
    Thus `filter_data_by_segment()` comes in handy, it will return a data frame that contains the rows that were logged during a time segment plus an extra column called `local_segment`. This new column will have as many unique values as time segment instances exist (14, 2, and 2 for our `p01`'s `my_days`, `my_weeks`, and `my_weekends` examples). After filtering, **you should group the data frame by this column and compute any desired features**, for example:

    ```python
    acc_features["maxmagnitude"] = acc_data.groupby(["local_segment"])["magnitude"].max()
    ```

    The reason RAPIDS does not filter the participant's data set for you is because your code might need to compute something based on a participant's complete dataset before computing their features. For example, you might want to identify the number that called a participant the most throughout the study before computing a feature with the number of calls the participant received from that number.

??? info "3. Return a data frame with your features"
    After filtering, grouping your data, and computing your features, your provider function should return a data frame that has:
    
    -  One row per time segment instance (e.g., 14 our `p01`'s `my_days` example)
    -  The `local_segment` column added by `filter_data_by_segment()`
    -  One column per feature. The name of your features should only contain letters or numbers (`feature1`) by convention. RAPIDS automatically adds the correct sensor and provider prefix; in our example, this prefix is `phone_accelerometr_vega_`.

??? example "`PHONE_ACCELEROMETER` Provider Example"
    For your reference, this our own provider (`RAPIDS`) for `PHONE_ACCELEROMETER` that computes five acceleration features

    ```python

    --8<---- "src/features/phone_accelerometer/rapids/main.py"

    ```

## New Features for Non-Existing Sensors

If you want to add features for a device or a sensor that we do not support at the moment (those that do not appear in the `"Existing Sensors"` list above), [open a new discussion](https://github.com/carissalow/rapids/discussions) in Github and we can add the necessary code so you can follow the instructions above.