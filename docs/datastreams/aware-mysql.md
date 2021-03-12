# `aware_mysql`

This [data stream](../../datastreams/data-streams-introduction) handles iOS and Android sensor data collected with the [AWARE Framework](https://awareframework.com/) and stored in a MySQL database.

## Container
A MySQL database with a table per sensor, each containing the data for all participants. This is the default database created by the old PHP AWARE server (as opposed to the new JavaScript Micro server).

The script to connect and download data from this container is at:
```bash
src/data/streams/aware_mysql/container.R
```

## Format

--8<---- "docs/snippets/aware_format.md"
