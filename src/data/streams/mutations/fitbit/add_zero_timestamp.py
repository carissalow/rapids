import pandas as pd

def main(parsed_data, stream_parameters):
    parsed_data["timestamp"] = 0 # this column is added at readable_datetime.R because we neeed to take into account multiple timezones
    # parsed_data['local_date_time'] = parsed_data['local_date_time'].dt.strftime('%Y-%m-%d %H:%M:%S')
    return(parsed_data)
