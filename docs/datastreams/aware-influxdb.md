# `aware_influxdb (beta)`

!!! warning
    This data stream is being released in beta while we test it thoroughly. 

This [data stream](../../datastreams/data-streams-introduction) handles iOS and Android sensor data collected with the [AWARE Framework](https://awareframework.com/) and stored in an InfluxDB database.

## Container
An InfluxDB database with a table per sensor, each containing the data for all participants.

The script to connect and download data from this container is at:
```bash
src/data/streams/aware_influxdb/container.R
```

## Format

--8<---- "docs/snippets/aware_format.md"
