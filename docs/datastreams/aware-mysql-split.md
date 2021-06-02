# `aware_mysql_split`

This [data stream](../../datastreams/data-streams-introduction) handles iOS and Android sensor data collected with the [AWARE Framework](https://awareframework.com/) and stored in a MySQL database. This stream is similar to `aware_mysql` except for the way data is stored in the database tables as explained below.

## Container
A MySQL database with a table per sensor **per participant**. RAPIDS assumes such tables' names follow the format `deviceid_sensorname` (for example `a748ee1a-1d0b-4ae9-9074-279a2b6ba524_accelerometer`); if this is not the case, you can modify the SQL query in this stream's `container.R`script. RAPIDS also assumes that an empty table exists for those participants that don’t have data for a specific sensor.

The script to connect and download data from this container is at:
```bash
src/data/streams/aware_mysql_split/container.R
```

## Format

--8<---- "docs/snippets/aware_format.md"
