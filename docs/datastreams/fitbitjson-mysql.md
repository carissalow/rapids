# `fitbitjson_mysql`
This [data stream](../../datastreams/data-streams-introduction) handles Fitbit sensor data downloaded using the [Fitbit Web API](https://dev.fitbit.com/build/reference/web-api/) and stored in a MySQL database. Please note that RAPIDS cannot query the API directly; you need to use other available tools or implement your own. Once you have your sensor data in a MySQL database, RAPIDS can process it.

## Container
The container should be a MySQL database with a table per sensor, each containing all participants' data.

The script to connect and download data from this container is at:
```bash
src/data/streams/fitbitjson_mysql/container.R
```

## Format

--8<---- "docs/snippets/jsonfitbit_format.md"
