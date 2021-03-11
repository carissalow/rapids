# `fitbitparsed_mysql`
This [data stream](../../datastreams/data-streams-introduction) handles Fitbit sensor data downloaded using the [Fitbit Web API](https://dev.fitbit.com/build/reference/web-api/), **parsed**, and stored in a MySQL database. Please note that RAPIDS cannot query the API directly; you need to use other available tools or implement your own. Once you have your parsed sensor data in a MySQL database, RAPIDS can process it.

!!! info "What is the difference between JSON and plain data streams"
    Most people will only need `fitbitjson_*` because they downloaded and stored their data directly from Fitbit's API. However, if, for some reason, you don't have access to that JSON data and instead only have the parsed data (columns and rows), you can use this data stream.

## Container
The container should be a MySQL database with a table per sensor, each containing all participants' data.

The script to connect and download data from this container is at:
```bash
src/data/streams/fitbitparsed_mysql/container.R
```

## Format

--8<---- "docs/snippets/parsedfitbit_format.md"
