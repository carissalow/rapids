# `aware_csv`

This [data stream](../../datastreams/data-streams-introduction) handles iOS and Android sensor data collected with the [AWARE Framework](https://awareframework.com/) and stored in CSV files.

!!! warning
    The CSV files have to use `,` as separator, `\` as escape character (do not escape `"` with `""`), and wrap any string columns with `"`.

    See examples in the CSV files inside [rapids_example_csv.zip](https://osf.io/wbg23/)

    ??? example "Example of a valid CSV file"
        ```csv
        "_id","timestamp","device_id","activities","confidence","stationary","walking","running","automotive","cycling","unknown","label"
        1,1587528000000,"13dbc8a3-dae3-4834-823a-4bc96a7d459d","[\"stationary\"]",2,1,0,0,0,0,0,""
        2,1587528060000,"13dbc8a3-dae3-4834-823a-4bc96a7d459d","[\"stationary\"]",2,1,0,0,0,0,0,"supplement"
        3,1587528120000,"13dbc8a3-dae3-4834-823a-4bc96a7d459d","[\"stationary\"]",2,1,0,0,0,0,0,"supplement"
        4,1587528180000,"13dbc8a3-dae3-4834-823a-4bc96a7d459d","[\"stationary\"]",2,1,0,0,0,0,0,"supplement"
        5,1587528240000,"13dbc8a3-dae3-4834-823a-4bc96a7d459d","[\"stationary\"]",2,1,0,0,0,0,0,"supplement"
        6,1587528300000,"13dbc8a3-dae3-4834-823a-4bc96a7d459d","[\"stationary\"]",2,1,0,0,0,0,0,"supplement"
        7,1587528360000,"13dbc8a3-dae3-4834-823a-4bc96a7d459d","[\"stationary\"]",2,1,0,0,0,0,0,"supplement"
        ```

## Container
A CSV file per sensor, each containing the data for all participants. 

The script to connect and download data from this container is at:
```bash
src/data/streams/aware_csv/container.R
```

## Format

--8<---- "docs/snippets/aware_format.md"
