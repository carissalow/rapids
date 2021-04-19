import pandas as pd

def main(parsed_data, stream_parameters):

    parsed_data["local_date_time"] = (pd.to_datetime(parsed_data["local_start_date_time"]) - pd.Timedelta(minutes=stream_parameters["SLEEP_SUMMARY_LAST_NIGHT_END"])).dt.strftime('%Y-%m-%d 00:00:00')

    # complete missing dates
    missed_dates = list(set([x.strftime('%Y-%m-%d 00:00:00') for x in pd.date_range(parsed_data["local_date_time"].min(), parsed_data["local_date_time"].max()).to_pydatetime()]) - set(parsed_data["local_date_time"]))
    parsed_data = pd.concat([parsed_data, pd.DataFrame({"local_date_time": missed_dates})], axis=0)
    parsed_data.sort_values(by=["local_date_time", "local_start_date_time"], inplace=True)
    parsed_data["device_id"] = parsed_data["device_id"].interpolate(method="pad")

    return parsed_data
