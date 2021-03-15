# `fitbitparsed_csv`
This [data stream](../../datastreams/data-streams-introduction) handles Fitbit sensor data downloaded using the [Fitbit Web API](https://dev.fitbit.com/build/reference/web-api/), **parsed**, and stored in a CSV file. Please note that RAPIDS cannot query the API directly; you need to use other available tools or implement your own. Once you have your parsed sensor data in a CSV file, RAPIDS can process it.

!!! info "What is the difference between JSON and plain data streams"
    Most people will only need `fitbitjson_*` because they downloaded and stored their data directly from Fitbit's API. However, if, for some reason, you don't have access to that JSON data and instead only have the parsed data (columns and rows), you can use this data stream.

!!! warning
    The CSV files have to use `,` as separator, `\` as escape character (do not escape `"` with `""`), and wrap any string columns with `"`.

    ??? example "Example of a valid CSV file"
        ```csv
        "device_id","heartrate","heartrate_zone","local_date_time","timestamp"
        "a748ee1a-1d0b-4ae9-9074-279a2b6ba524",69,"outofrange","2020-04-23 00:00:00",0
        "a748ee1a-1d0b-4ae9-9074-279a2b6ba524",69,"outofrange","2020-04-23 00:01:00",0
        "a748ee1a-1d0b-4ae9-9074-279a2b6ba524",67,"outofrange","2020-04-23 00:02:00",0
        "a748ee1a-1d0b-4ae9-9074-279a2b6ba524",69,"outofrange","2020-04-23 00:03:00",0
        ```

## Container
The container should be a CSV file per sensor, each containing all participants' data.

The script to connect and download data from this container is at:
```bash
src/data/streams/fitbitparsed_csv/container.R
```

## Format

--8<---- "docs/snippets/parsedfitbit_format.md"
