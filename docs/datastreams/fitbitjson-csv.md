# `fitbitjson_csv`
This [data stream](../../datastreams/data-streams-introduction) handles Fitbit sensor data downloaded using the [Fitbit Web API](https://dev.fitbit.com/build/reference/web-api/) and stored in a CSV file. Please note that RAPIDS cannot query the API directly; you need to use other available tools or implement your own. Once you have your sensor data in a CSV file, RAPIDS can process it.

## Container
The container should be a CSV file per Fitbit sensor, each containing all participants' data.

The script to connect and download data from this container is at:
```bash
src/data/streams/fitbitjson_csv/container.R
```

## Format

--8<---- "docs/snippets/jsonfitbit_format.md"
