import pandas as pd

def main(parsed_data, stream_parameters):

    if parsed_data.empty:
        return pd.DataFrame(columns=parsed_data.columns.tolist() + ["local_date_time"])
    
    parsed_data["local_date_time"] = (pd.to_datetime(parsed_data["local_start_date_time"]) - pd.Timedelta(minutes=stream_parameters["SLEEP_SUMMARY_LAST_NIGHT_END"]))
    start_date, end_date = parsed_data["local_date_time"].min().strftime('%Y-%m-%d 00:00:00'), (parsed_data["local_date_time"].max() + pd.Timedelta(days=1)).strftime('%Y-%m-%d 00:00:00')
    parsed_data["local_date_time"] = parsed_data["local_date_time"].dt.strftime('%Y-%m-%d 00:00:00')

    # complete missing dates
    missed_dates = list(set([x.strftime('%Y-%m-%d 00:00:00') for x in pd.date_range(start_date, end_date).to_pydatetime()]) - set(parsed_data["local_date_time"]))
    parsed_data = pd.concat([parsed_data, pd.DataFrame({"local_date_time": missed_dates})], axis=0)
    parsed_data.sort_values(by=["local_date_time", "local_start_date_time"], inplace=True)
    parsed_data["device_id"] = parsed_data["device_id"].interpolate(method="pad")

    if pd.api.types.is_datetime64_any_dtype(parsed_data['local_start_date_time']):
        parsed_data['local_start_date_time'] = parsed_data['local_start_date_time'].dt.strftime('%Y-%m-%d %H:%M:%S')
    if pd.api.types.is_datetime64_any_dtype(parsed_data['local_end_date_time']):
        parsed_data['local_end_date_time'] = parsed_data['local_end_date_time'].dt.strftime('%Y-%m-%d %H:%M:%S')

    return parsed_data
