# `aware_micro_mysql`

This [data stream](../../datastreams/data-streams-introduction) handles iOS and Android sensor data collected with the [AWARE Framework's](https://awareframework.com/) [AWARE Micro](https://github.com/denzilferreira/aware-micro) server and stored in a MySQL database.

## Container
A MySQL database with a table per sensor, each containing the data for all participants. Sensor data is stored in a JSON field within each table called `data`

The script to connect and download data from this container is at:
```bash
src/data/streams/aware_micro_mysql/container.R
```

## Format

--8<---- "docs/snippets/aware_format.md"