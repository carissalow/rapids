# Phone Log

Sensor parameters description for `[PHONE_LOG]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER][ANDROID]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where a data log is stored for Android devices
|`[CONTAINER][IOS]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where a data log is stored for iOS devices

!!! note
    No feature providers have been implemented for this sensor yet, however you can use its key (`PHONE_LOG`) to improve [`PHONE_DATA_YIELD`](../phone-data-yield) or you can [implement your own features](../add-new-features).