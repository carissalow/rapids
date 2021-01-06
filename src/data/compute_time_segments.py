import pandas as pd
import warnings
import yaml

def is_valid_frequency_segments(time_segments, time_segments_file):
    """
    returns true if time_segment has the expected structure for generating frequency segments;
    raises ValueError exception otherwise.
    """
    
    valid_columns = ["label", "length"]
    if set(time_segments.columns) != set(valid_columns):
        error_message = 'The FREQUENCY time segments file in [TIME_SEGMENTS][FILE] must have two columns: label, and length ' \
                  'but instead we found {}. Modify {}'.format(list(time_segments.columns), time_segments_file)
        raise ValueError(error_message)

    if time_segments.shape[0] > 1:
        message = 'The FREQUENCY time segments file in [TIME_SEGMENTS][FILE] can only have 1 row.' \
                  'Modify {}'.format(time_segments_file)
        raise ValueError(message)

    if not pd.api.types.is_integer_dtype(time_segments.dtypes['length']):
        message = 'The column length in the FREQUENCY time segments file in [TIME_SEGMENTS][FILE] must be integer but instead is ' \
                  '{}. . This usually means that not all values in this column are formed by digits. Modify {}'.format(time_segments.dtypes['length'], time_segments_file)
        raise ValueError(message)

    if time_segments.iloc[0].loc['length'] < 0:
        message = 'The value in column length in the FREQUENCY time segments file in [TIME_SEGMENTS][FILE] must be positive but instead is  ' \
                  '{}. Modify {}'.format(time_segments.iloc[0].loc['length'], time_segments_file)
        raise ValueError(message)
    if time_segments.iloc[0].loc['length'] >= 1440:
        message = 'The column length in the FREQUENCY time segments file in [TIME_SEGMENTS][FILE] must be shorter than a day in minutes (1440) but instead is ' \
                  '{}. Modify {}'.format(time_segments.iloc[0].loc['length'], time_segments_file)
        raise ValueError(message)

    return True

def is_valid_periodic_segments(time_segments, time_segments_file):
    time_segments = time_segments.copy(deep=True)

    valid_columns = ["label", "start_time", "length", "repeats_on", "repeats_value"]
    if set(time_segments.columns) != set(valid_columns):
        error_message = 'The PERIODIC time segments file in [TIME_SEGMENTS][FILE] must have five columns: label, start_time, length, repeats_on, repeats_value ' \
                  'but instead we found {}. Modify {}'.format(list(time_segments.columns), time_segments_file)
        raise ValueError(error_message)

    valid_repeats_on = ["every_day", "wday", "mday", "qday", "yday"]
    if len(list(set(time_segments["repeats_on"]) - set(valid_repeats_on))) > 0:
        error_message = 'The column repeats_on in the PERIODIC time segments file in [TIME_SEGMENTS][FILE] can only accept: "every_day", "wday", "mday", "qday", or "yday" ' \
                  'but instead we found {}. Modify {}'.format(list(set(time_segments["repeats_on"])), time_segments_file)
        raise ValueError(error_message)

    if not pd.api.types.is_integer_dtype(time_segments.dtypes['repeats_value']):
        message = 'The column repeats_value in the PERIODIC time segments file in [TIME_SEGMENTS][FILE] must be integer but instead is ' \
                  '{}. . This usually means that not all values in this column are formed by digits. Modify {}'.format(time_segments.dtypes['repeats_value'], time_segments_file)
        raise ValueError(message)

    invalid_time_segments = time_segments.query("repeats_on == 'every_day' and repeats_value != 0")
    if invalid_time_segments.shape[0] > 0:
        message = 'Every row with repeats_on=every_day must have a repeats_value=0 in the PERIODIC time segments file in [TIME_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_time_segments["label"].to_numpy(), time_segments_file)
        raise ValueError(message)

    invalid_time_segments = time_segments.query("repeats_on == 'wday' and (repeats_value < 1 | repeats_value > 7)")
    if invalid_time_segments.shape[0] > 0:
        message = 'Every row with repeats_on=wday must have a repeats_value=[1,7] in the PERIODIC time segments file in [TIME_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_time_segments["label"].to_numpy(), time_segments_file)
        raise ValueError(message)

    invalid_time_segments = time_segments.query("repeats_on == 'mday' and (repeats_value < 1 | repeats_value > 31)")
    if invalid_time_segments.shape[0] > 0:
        message = 'Every row with repeats_on=mday must have a repeats_value=[1,31] in the PERIODIC time segments file in [TIME_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_time_segments["label"].to_numpy(), time_segments_file)
        raise ValueError(message)

    invalid_time_segments = time_segments.query("repeats_on == 'qday' and (repeats_value < 1 | repeats_value > 92)")
    if invalid_time_segments.shape[0] > 0:
        message = 'Every row with repeats_on=qday must have a repeats_value=[1,92] in the PERIODIC time segments file in [TIME_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_time_segments["label"].to_numpy(), time_segments_file)
        raise ValueError(message)

    invalid_time_segments = time_segments.query("repeats_on == 'yday' and (repeats_value < 1 | repeats_value > 366)")
    if invalid_time_segments.shape[0] > 0:
        message = 'Every row with repeats_on=yday must have a repeats_value=[1,366] in the PERIODIC time segments file in [TIME_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_time_segments["label"].to_numpy(), time_segments_file)
        raise ValueError(message)

    try:
        time_segments["start_time"] = pd.to_datetime(time_segments["start_time"])
    except ValueError as err:
        raise ValueError("At least one start_time in the PERIODIC time segments file in [TIME_SEGMENTS][FILE] has an invalid format, it should be HH:MM:SS in 24hr clock({}). Modify {}".format(err, time_segments_file))

    if(time_segments.shape[0] != time_segments.drop_duplicates().shape[0]):
        error_message = 'The PERIODIC time segments file in [TIME_SEGMENTS][FILE] has two or more rows that are identical. ' \
                  'Modify {}'.format(time_segments_file)
        raise ValueError(error_message)
    
    duplicated_labels = time_segments[time_segments["label"].duplicated()]
    if(duplicated_labels.shape[0] > 0):
        error_message = 'Segements labels must be unique. The PERIODIC time segments file in [TIME_SEGMENTS][FILE] has {} row(s) with the same label {}. ' \
                  'Modify {}'.format(duplicated_labels.shape[0], duplicated_labels["label"].to_numpy(), time_segments_file)
        raise ValueError(error_message)
    
    # TODO Validate string format for lubridate
    
    return True

def is_valid_event_segments(time_segments, time_segments_file):
    time_segments = time_segments.copy(deep=True)

    valid_columns = ["label", "event_timestamp", "length", "shift", "shift_direction", "device_id"]
    if set(time_segments.columns) != set(valid_columns):
        error_message = 'The EVENT time segments file in [TIME_SEGMENTS][FILE] must have six columns: label, event_timestamp, length, shift, shift_direction and device_id ' \
                  'but instead we found {}. Modify {}'.format(list(time_segments.columns), time_segments_file)
        raise ValueError(error_message)

    if not pd.api.types.is_integer_dtype(time_segments.dtypes['event_timestamp']):
        message = 'The column event_timestamp in the EVENT time segments file in [TIME_SEGMENTS][FILE] must be integer but instead is ' \
                  '{}. This usually means that not all values in this column are formed by digits. Modify {}'.format(time_segments.dtypes['event_timestamp'], time_segments_file)
        raise ValueError(message)

    valid_shift_direction_values = [1, -1, 0]
    provided_values = time_segments["shift_direction"].unique()
    if len(list(set(provided_values) - set(valid_shift_direction_values))) > 0:
        error_message = 'The values of shift_direction column in the EVENT time segments file in [TIME_SEGMENTS][FILE] can only be 1, -1 or 0 ' \
                  'but instead we found {}. Modify {}'.format(provided_values, time_segments_file)
        raise ValueError(error_message)

    if(time_segments.shape[0] != time_segments.drop_duplicates().shape[0]):
        error_message = 'The EVENT time segments file in [TIME_SEGMENTS][FILE] has two or more rows that are identical. ' \
                  'Modify {}'.format(time_segments_file)
        raise ValueError(error_message)

    # TODO Validate string format for lubridate of length and shift
    # TODO validate unique labels per participant
    return True


def parse_frequency_segments(time_segments: pd.DataFrame) -> pd.DataFrame:
    """
    returns a table with rows identifying start and end of time slots with frequency freq (in minutes). For example,
    for freq = 10 it outputs:
        bin_id start end   label
        0      00:00 00:10 epoch_0000
        1      00:10 00:20 epoch_0001
        2      00:20 00:30 epoch_0002
        ...
        143    23:50 00:00 epoch_0143
    time_segments argument is expected to have the following structure:
        label  length
        epoch      10
    """
    freq = time_segments.iloc[0].loc['length']
    slots = pd.date_range(start='2020-01-01', end='2020-01-02', freq='{}min'.format(freq))
    slots = ['{:02d}:{:02d}'.format(x.hour, x.minute) for x in slots]

    table = pd.DataFrame(slots, columns=['start_time'])
    table['length'] = time_segments.iloc[0].loc['length']
    table = table.iloc[:-1, :]

    label = time_segments.loc[0, 'label']
    table['label'] = range(0, table.shape[0])
    table['label'] = table['label'].apply(lambda x: '{}{:04}'.format(label, x))

    return table[['start_time', 'length', 'label']]

def parse_periodic_segments(time_segments):
    time_segments.loc[time_segments["repeats_on"] == "every_day", "repeats_value"] = 0
    return time_segments

def parse_event_segments(time_segments, device_ids):
    return time_segments.query("device_id == @device_ids")

def parse_time_segments(time_segments_file, segments_type, device_ids):
    # Add code to validate and parse frequencies, intervals, and events
    # Expected formats:
    # Frequency: label, length columns (e.g. my_prefix, 5) length has to be in minutes (int)
    # Interval: label, start, end columns (e.g. daily, 00:00, 23:59) start and end should be valid hours in 24 hour format
    # Event: label, timestamp, length, shift (e.g., survey1, 1532313215463, 60, -30), timestamp is a UNIX timestamp in ms (we could take a date time string instead), length is in minutes (int), shift is in minutes (+/-int) and is added/substracted from timestamp
    # Our output should have local_date, start_time, end_time, label. In the readable_datetime script, If local_date has the same value for all rows, every segment will be applied for all days, otherwise each segment will be applied only to its local_date
    time_segments = pd.read_csv(time_segments_file)

    if time_segments is None:
        message = 'The time segments file in [TIME_SEGMENTS][FILE] is None. Modify {}'.format(time_segments_file)
        raise ValueError(message)

    if time_segments.shape[0] == 0:
        message = 'The time segments file in [TIME_SEGMENTS][FILE] is empty. Modify {}'.format(time_segments_file)
        raise ValueError(message)

    if(segments_type not in ["FREQUENCY", "PERIODIC", "EVENT"]):
        raise ValueError("[TIME_SEGMENTS][TYPE] can only be FREQUENCY, PERIODIC, or EVENT")
    
    if(segments_type == "FREQUENCY" and is_valid_frequency_segments(time_segments, time_segments_file)):
        time_segments = parse_frequency_segments(time_segments)
    elif(segments_type == "PERIODIC" and is_valid_periodic_segments(time_segments, time_segments_file)):
        time_segments = parse_periodic_segments(time_segments)
    elif(segments_type == "EVENT" and is_valid_event_segments(time_segments, time_segments_file)):
        time_segments = parse_event_segments(time_segments, device_ids)
    else:
        raise ValueError("{} does not have a format compatible with frequency, periodic or event time segments. Please refer to [LINK]".format(time_segments_file))
    return time_segments

participant_file = yaml.load(open(snakemake.input[1], 'r'), Loader=yaml.FullLoader)
device_ids = []
for key in participant_file.keys():
    if "DEVICE_IDS" in participant_file[key] and isinstance(participant_file[key]["DEVICE_IDS"], list):
        device_ids = device_ids + participant_file[key]["DEVICE_IDS"]

final_time_segments = parse_time_segments(snakemake.input[0], snakemake.params["time_segments_type"], device_ids)

if snakemake.params["time_segments_type"] == "EVENT" and final_time_segments.shape[0] == 0:
    warnings.warn("There are no event time segments for {}. Check your time segment file {}".format(snakemake.params["pid"], snakemake.input[0]))

final_time_segments.to_csv(snakemake.output["segments_file"], index=False)
pd.DataFrame({"label" : final_time_segments["label"].unique()}).to_csv(snakemake.output["segments_labels_file"], index=False)