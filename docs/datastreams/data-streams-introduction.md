# Data Streams Introduction

A data stream is a set of sensor data collected using a specific type of **device** with a specific **format** and stored in a specific **container**.

For example, the `aware_mysql` data stream handles smartphone data (**device**) collected with the [AWARE Framework](https://awareframework.com/) (**format**) stored in a MySQL database (**container**). Similarly, smartphone data collected with [Beiwe](https://www.beiwe.org/) will have a different format and could be stored in a container like a PostgreSQL database or a CSV file.

If you want to process a data stream using RAPIDS, make sure that your data is stored in a supported **format** and **container** (see table below). 

If RAPIDS doesn't support your data stream yet (e.g. Beiwe data stored in PostgreSQL, or AWARE data stored in SQLite), you can always [implement a new data stream](../add-new-data-streams). If it's something you think other people might be interested on, we will be happy to include your new data stream in RAPIDS, so get in touch!.

!!! hint
    Currently, you can add new data streams for smartphones, Fitbit, and Empatica devices. If you need RAPIDS to process data from **other devices**, like Oura Rings or Actigraph wearables, get in touch. It is a more complicated process that could take a couple of days to implement for someone familiar with R or Python, but we would be happy to work on it together.

For reference, these are the data streams we currently support: 

| Data Stream | Device | Format | Container | Docs
|--|--|--|--|--|
| `aware_mysql`| Phone | AWARE app | MySQL | [link](../aware-mysql)
| `aware_micro_mysql`| Phone | AWARE Micro server | MySQL | [link](../aware-micro-mysql)
| `aware_csv`| Phone | AWARE app | CSV files | [link](../aware-csv)
| `aware_influxdb` (beta)| Phone | AWARE app | InfluxDB | [link](../aware-influxdb)
| `fitbitjson_mysql`| Fitbit | JSON (per [Fitbit's API](https://dev.fitbit.com/build/reference/web-api/)) | MySQL | [link](../fitbitjson-mysql)
| `fitbitjson_csv`| Fitbit | JSON (per [Fitbit's API](https://dev.fitbit.com/build/reference/web-api/)) | CSV files | [link](../fitbitjson-csv)
| `fitbitparsed_mysql`| Fitbit | Parsed (parsed API data) | MySQL | [link](../fitbitparsed-mysql)
| `fitbitparsed_csv`| Fitbit | Parsed (parsed API data)  | CSV files | [link](../fitbitparsed-csv)
| `empatica_zip`| Empatica | [E4 Connect](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-) | ZIP files | [link](../empatica-zip)
