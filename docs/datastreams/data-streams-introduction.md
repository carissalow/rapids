# Data Streams Introduction

A data stream is a set of sensor data collected using a specific type of **device** with a specific **format** and stored in a specific **container**.

For example, the `aware_mysql` data stream handles smartphone data (**device**) collected with the [AWARE Framework](https://awareframework.com/) (**format**) stored in a MySQL database (**container**). Similarly, smartphone data collected with [Beiwe](https://www.beiwe.org/) will have a different format and could be stored in a container like a PostgreSQL database or a CSV file.

If you want to process a data stream using RAPIDS, make sure that your data is stored in a supported **format** and **container** (see table below). 

If RAPIDS doesn't support your data stream yet (e.g. Beiwe data stored in PostgreSQL, or AWARE data stored in InfluxDB), you can always [implement a new data stream](../add-new-data-streams). If it's something you think other people might be interested on, we will be happy to include your new data stream in RAPIDS, so get in touch!.

!!! hint
    You can only add new data streams for Smartphone or Fitbit data. If you need RAPIDS to process data from **different devices**, like Oura Rings or Actigraph wearables, get in touch. It is a more complex process that could take a few days to implement for someone familiar with R or Python but that we would be happy to work on together.

For reference, these are the data streams we currently support: 

| Data Stream | Device | Format | Container | Docs
|--|--|--|--|--|
| `aware_mysql`| Phone | AWARE app | MySQL | [link]()
| `aware_csv`| Phone | AWARE app | CSV files | [link]()
| `fitbitjson_mysql`| Fitbit | JSON (per Fitbit's API) | MySQL | [link]()
| `fitbitjson_csv`| Fitbit | JSON (per Fitbit's API) | CSV files | [link]()
| `fitbitparsed_mysql`| Fitbit | Parsed (parsed API data) | MySQL | [link]()
| `fitbitparsed_csv`| Fitbit | Parsed (parsed API data)  | CSV files | [link]()
| `empatica_zip`| Empatica | E4 Connect | ZIP files | [link]()

!!! hint
    - Fitbit data can be processed from the JSON object produced by Fitbit's API (recommended) or from parsed tabular data (if you only have access to parsed data).
    - Empatica data can only be accessed through the [E4 Connect website](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-) that produces zip files with a CSV file per sensor which can be processed directly in RAPIDS. 