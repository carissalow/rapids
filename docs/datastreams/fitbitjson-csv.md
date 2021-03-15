# `fitbitjson_csv`
This [data stream](../../datastreams/data-streams-introduction) handles Fitbit sensor data downloaded using the [Fitbit Web API](https://dev.fitbit.com/build/reference/web-api/) and stored in a CSV file. Please note that RAPIDS cannot query the API directly; you need to use other available tools or implement your own. Once you have your sensor data in a CSV file, RAPIDS can process it.

!!! warning
    The CSV files have to use `,` as separator, `\` as escape character (do not escape `"` with `""`), and wrap any string columns with `"`.

    ??? example "Example of a valid CSV file"
        ```csv
        "timestamp","device_id","label","fitbit_id","fitbit_data_type","fitbit_data"
        1587614400000,"a748ee1a-1d0b-4ae9-9074-279a2b6ba524","5S","5ZKN9B","steps","{\"activities-steps\":[{\"dateTime\":\"2020-04-23\",\"value\":\"7881\"}]"
        ```

## Container
The container should be a CSV file per Fitbit sensor, each containing all participants' data.

The script to connect and download data from this container is at:
```bash
src/data/streams/fitbitjson_csv/container.R
```

## Format

--8<---- "docs/snippets/jsonfitbit_format.md"
