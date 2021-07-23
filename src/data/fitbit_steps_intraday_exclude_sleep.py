import pandas as pd
from datetime import timedelta


exclude_sleep = snakemake.params["exclude_sleep"]
exclude_time_based = exclude_sleep["TIME_BASED"]["EXCLUDE"]
exclude_fitbit_based = exclude_sleep["FITBIT_BASED"]["EXCLUDE"]
exclude_sleep_fixed_start = exclude_sleep["TIME_BASED"]["START_TIME"] + ":00"
exclude_sleep_fixed_end = exclude_sleep["TIME_BASED"]["END_TIME"] + ":00"


steps_intraday_data = pd.read_csv(snakemake.input["sensor_data"], parse_dates=["local_date"])
sleep_data = pd.read_csv(snakemake.input["sleep_data"]).dropna(subset=["local_start_date_time", "local_end_date_time"], how="any") if snakemake.input["sleep_data"] else pd.DataFrame()

if not steps_intraday_data.empty:

    if exclude_time_based and (not exclude_fitbit_based or (exclude_fitbit_based and sleep_data.empty)):
        query_string = "local_time < @exclude_sleep_fixed_start " + ("&" if exclude_sleep_fixed_start >= exclude_sleep_fixed_end else "|") + " local_time > @exclude_sleep_fixed_end"
        steps_intraday_data.query(query_string, inplace=True)

    elif exclude_fitbit_based and (not sleep_data.empty):

        queries = []

        if exclude_time_based:

            # Get fixed intervals
            fixed_start_dates = pd.date_range(steps_intraday_data["local_date"].min() - timedelta(days=1), steps_intraday_data["local_date"].max())
            fixed_end_dates = fixed_start_dates + timedelta(days=1) if exclude_sleep_fixed_start >= exclude_sleep_fixed_end else fixed_start_dates

            fixed_time = pd.DataFrame({"fixed_start_date_time": fixed_start_dates.strftime("%Y-%m-%d") + " " + exclude_sleep_fixed_start,
                            "fixed_end_date_time": fixed_end_dates.strftime("%Y-%m-%d") + " " + exclude_sleep_fixed_end})
            
            # Remove fixed intervals that intersect with sleep intervals from the fixed sleep periods
            sleep_data["query_intersect"] = "(fixed_start_date_time < '" + sleep_data["local_end_date_time"] + "' & fixed_end_date_time > '" + sleep_data["local_start_date_time"] + "')"
            query_string_subtract = "~(" + " | ".join(sleep_data["query_intersect"].tolist()) + ")"
            fixed_time.query(query_string_subtract, inplace=True)
            
            # Add TIME_BASED query to queries
            queries = ("(local_date_time < '" + fixed_time["fixed_start_date_time"] + "' | local_date_time > '" + fixed_time["fixed_end_date_time"] + "')").tolist()

        # Add FITBIT_BASED query to queries
        queries = queries + ("(local_date_time < '" + sleep_data["local_start_date_time"] + "' | local_date_time > '" + sleep_data["local_end_date_time"] + "')").tolist()

        steps_intraday_data.query(" & ".join(queries), inplace=True)

steps_intraday_data.to_csv(snakemake.output[0], index=False)
