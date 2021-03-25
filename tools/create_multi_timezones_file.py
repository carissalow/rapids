"""
This script can create the multiple timezones file based on timezone table collected by the AWARE app.

Input: timezone table collected by the AWARE app
---
Expected output:

| Column      | Description                                                                                                   |
|-------------|---------------------------------------------------------------------------------------------------------------|
| device_id   | A string that uniquely identifies a smartphone or wearable                                                    |
| tzcode      | A string with the appropriate code from this list that represents the time zone where the device sensed data  |
| timestamp   | A UNIX timestamp indicating when was the first time this device_id sensed data in tzcode                      |

How to run it?
1. Put the timezone table (timezone.csv) collected by the AWARE app under data/external folder
2. Run python tools/create_multi_timezones_file.py

"""

import re
import pandas as pd

# Load the timezone table collected by the AWARE app
data = pd.read_csv("data/external/timezone.csv")
# Load the first table of the List of tz database time zones page in wiki (Link: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones)
wiki_tz = pd.read_csv("data/external/wiki_tz.csv")

data["tzcode"] = data["timezone"].str.extract("("+"|".join(wiki_tz["TZ database name"].tolist())+")", expand=False)
data = data[["device_id", "tzcode", "timestamp"]]

# Sort by device_id and timestamp
data.sort_values(by=["device_id", "timestamp"], inplace=True)
# Only keep the first & last row for consecutive rows with the same tzcode per device_id
data_first = data.loc[(data["device_id"].shift(1) != data["device_id"]) | (data["tzcode"].shift(1) != data["tzcode"])]
data_last = data.loc[(data["device_id"].shift(-1) != data["device_id"]) | (data["tzcode"].shift(-1) != data["tzcode"])]
data = pd.concat([data_first, data_last], axis=0)
# Drop duplicates and sort by device_id and timestamp
data = data.drop_duplicates()
data.sort_values(by=["device_id", "timestamp"], inplace=True)

data.to_csv("data/external/multiple_timezones.csv", index=False)
