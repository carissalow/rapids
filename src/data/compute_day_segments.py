import pandas as pd
import json

def parse_day_segments(segments, event_time_shift, event_segment_duration):
    # Temporal code to parse segments, should substitute with the code to parse
    # frequencies, intervals, and events
    data = json.loads(segments)
    label = []
    start = []
    end = []
    for d in data:
        label.append(d[0])
        start.append(d[1])
        end.append(d[2])

    day_segments = pd.DataFrame(list(zip([1]*len(label), start, end, label)), columns =['local_date','start_time','end_time','label']) 
    return day_segments
    ##########################

segments = snakemake.params["segments"]
event_time_shift = snakemake.params["event_time_shift"]
event_segment_duration = snakemake.params["event_segment_duration"]

day_segments = parse_day_segments(segments, event_time_shift, event_segment_duration)
day_segments.to_csv(snakemake.output[0], index=False)