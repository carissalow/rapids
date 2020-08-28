
def filter_data_by_segment(data, day_segment):
    date_regex = "[0-9]{4}[\-|\/][0-9]{2}[\-|\/][0-9]{2}"
    hour_regex = "[0-9]{2}:[0-9]{2}:[0-9]{2}"
    segment_regex = "\[({}#{}#{}#{}#{})\]".format(day_segment, date_regex, hour_regex, date_regex, hour_regex)
    data["local_segment"] = data["assigned_segments"].str.extract(segment_regex, expand=True)
    return(data.dropna(subset = ["local_segment"]))

rapids_log_tag =  "RAPIDS:"