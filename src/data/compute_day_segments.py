import pandas as pd


def parse_day_segments(day_segments):
    # Add code to parse frequencies, intervals, and events
    # Expected formats:
    # Frequency: label, length columns (e.g. my_prefix, 5) length has to be in minutes (int)
    # Interval: label, start, end columns (e.g. daily, 00:00, 23:59) start and end should be valid hours in 24 hour format
    # Event: label, timestamp, length, shift (e.g., survey1, 1532313215463, 60, -30), timestamp is a UNIX timestamp in ms, length is in minutes (int), shift is in minutes (+/-int) and is add/substracted from timestamp
    # Our output should have local_date, start_time, end_time, label. In the readable_datetime script, If local_date has the same value for all rows, every segment will be applied for all days, otherwise each segment will be applied only to its local_date
    
    day_segments["local_date"] = 1
    day_segments = day_segments.rename(columns={"start": "start_time", "end":"end_time"})
    return day_segments
    ##########################

day_segments = pd.read_csv(snakemake.input[0])
day_segments = parse_day_segments(day_segments)
day_segments.to_csv(snakemake.output["segments_file"], index=False)