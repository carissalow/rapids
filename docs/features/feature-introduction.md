# Behavioral Features Introduction

Every phone or Fitbit sensor has a corresponding config section e.g. `[PHONE_CALLS]` in `config.yaml` and 0, 1 or more **providers**. A provider is a script created by the core RAPIDS team or other researchers that extracts behavioral features for that sensor.

If you want to extract features from any sensor, set the corresponding `[PROVIDER][COMPUTE]` variable to `TRUE`, the `[TABLE]` variable to the sensor's table name in your database, and change any other parameters as [desired](/setup/configuration/#sensor-and-features-to-process), and [execute](/setup/execution/) RAPIDS:

!!! example
    In this example the `config.yaml` file has been configured to extract `PHONE_CALLS` features from a table called `calls`:
    ```yaml
    PHONE_CALLS:
        TABLE: calls
        PROVIDERS:
            RAPIDS:
                COMPUTE: True
            ...
    ```

!!! hint
    Every time you change any sensor parameters, all the necessary files will be updated as soon as you execute RAPIDS. Some sensors will have specific attributes (like `MESSAGES_TYPES`) so refer to each sensor documentation.

Click on the left-hand side menu for the full feature catalogue per sensor.
