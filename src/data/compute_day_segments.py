import pandas as pd

def is_valid_frequency_segments(day_segments, day_segments_file):
    """
    returns true if day_segment has the expected structure for generating frequency segments;
    raises ValueError exception otherwise.
    """
    
    valid_columns = ["label", "length"]
    if len(list(set(day_segments.columns) - set(valid_columns))) > 0:
        error_message = 'The FREQUENCY_EVERY_DAY day segments file in [DAY_SEGMENTS][FILE] must have two columns: label, and length ' \
                  'but instead we found {}. Modify {}'.format(list(day_segments.columns), day_segments_file)
        raise ValueError(error_message)

    if day_segments.shape[0] > 1:
        message = 'The FREQUENCY_EVERY_DAY day segments file in [DAY_SEGMENTS][FILE] can only have 1 row.' \
                  'Modify {}'.format(day_segments_file)
        raise ValueError(message)

    if not pd.api.types.is_integer_dtype(day_segments.dtypes['length']):
        message = 'The column length in the FREQUENCY_EVERY_DAY day segments file in [DAY_SEGMENTS][FILE] must be integer but instead is ' \
                  '{}. Modify {}'.format(day_segments.dtypes['length'], day_segments_file)
        raise ValueError(message)

    if day_segments.iloc[0].loc['length'] < 0:
        message = 'The value in column length in the FREQUENCY_EVERY_DAY day segments file in [DAY_SEGMENTS][FILE] must be positive but instead is  ' \
                  '{}. Modify {}'.format(day_segments.iloc[0].loc['length'], day_segments_file)
        raise ValueError(message)
    if day_segments.iloc[0].loc['length'] >= 1440:
        message = 'The column length in the FREQUENCY_EVERY_DAY day segments file in [DAY_SEGMENTS][FILE] must be shorter than a day in minutes (1440) but instead is ' \
                  '{}. Modify {}'.format(day_segments.iloc[0].loc['length'], day_segments_file)
        raise ValueError(message)

    return True

def is_valid_interval_segments(day_segments, day_segments_file):
    day_segments = day_segments.copy(deep=True)

    valid_columns = ["label", "start_time", "length"]
    if len(list(set(day_segments.columns) - set(valid_columns))) > 0:
        error_message = 'The INTERVAL_EVERY_DAY day segments file in [DAY_SEGMENTS][FILE] must have three columns: label, start_time and length ' \
                  'but instead we found {}. Modify {}'.format(list(day_segments.columns), day_segments_file)
        raise ValueError(error_message)

    try:
        day_segments["start_time"] = pd.to_datetime(day_segments["start_time"])
    except ValueError as err:
        raise ValueError("At least one start_time in the INTERVAL_EVERY_DAY day segments file in [DAY_SEGMENTS][FILE] has an invalid format, it should be HH:MM in 24hr clock({}). Modify {}".format(err, day_segments_file))

    if(day_segments.shape[0] != day_segments.drop_duplicates().shape[0]):
        error_message = 'The INTERVAL_EVERY_DAY day segments file in [DAY_SEGMENTS][FILE] has two or more rows that are identical. ' \
                  'Modify {}'.format(day_segments_file)
        raise ValueError(error_message)

    # TODO Validate string format for lubridate
    
    return True

def is_valid_event_segments(day_segments, day_segments_file):
    day_segments = day_segments.copy(deep=True)

    valid_columns = ["label", "start_date_time", "length", "shift", "shift_direction"]
    if len(list(set(day_segments.columns) - set(valid_columns))) > 0:
        error_message = 'The INTERVAL_FLEXIBLE_DAY day segments file in [DAY_SEGMENTS][FILE] must have five columns: label, start_date_time, length, shift and shift_direction ' \
                  'but instead we found {}. Modify {}'.format(list(day_segments.columns), day_segments_file)
        raise ValueError(error_message)

    try:
        day_segments["start_date_time"] = pd.to_datetime(day_segments["start_date_time"], format='%Y-%m-%d %H:%M:%S', errors='raise')
    except ValueError as err:
        raise ValueError("At least one start_date_time has an invalid format, it should be YYYY-MM-DD HH:MM:SS in 24hr clock({}). Modify {}".format(err, day_segments_file))

    valid_shift_direction_values = [1, -1, 0]
    provided_values = day_segments["shift_direction"].unique()
    if len(list(set(provided_values) - set(valid_shift_direction_values))) > 0:
        error_message = 'The values of shift_direction column in the INTERVAL_FLEXIBLE_DAY day segments file in [DAY_SEGMENTS][FILE] can only be 1, -1 or 0 ' \
                  'but instead we found {}. Modify {}'.format(provided_values, day_segments_file)
        raise ValueError(error_message)

    if(day_segments.shape[0] != day_segments.drop_duplicates().shape[0]):
        error_message = 'The INTERVAL_FLEXIBLE_DAY day segments file in [DAY_SEGMENTS][FILE] has two or more rows that are identical. ' \
                  'Modify {}'.format(day_segments_file)
        raise ValueError(error_message)

    # TODO Validate string format for lubridate of length and shift
    return True


def parse_frequency_segments(day_segments: pd.DataFrame) -> pd.DataFrame:
    """
    returns a table with rows identifying start and end of time slots with frequency freq (in minutes). For example,
    for freq = 10 it outputs:
        bin_id start end   label
        0      00:00 00:10 epoch_0000
        1      00:10 00:20 epoch_0001
        2      00:20 00:30 epoch_0002
        ...
        143    23:50 00:00 epoch_0143
    day_segments argument is expected to have the following structure:
        label  length
        epoch      10
    """
    freq = day_segments.iloc[0].loc['length']
    slots = pd.date_range(start='2020-01-01', end='2020-01-02', freq='{}min'.format(freq))
    slots = ['{:02d}:{:02d}'.format(x.hour, x.minute) for x in slots]

    table = pd.DataFrame(slots, columns=['start_time'])
    table['length'] = day_segments.iloc[0].loc['length']
    table = table.iloc[:-1, :]

    label = day_segments.loc[0, 'label']
    table['label'] = range(0, table.shape[0])
    table['label'] = table['label'].apply(lambda x: '{}{:04}'.format(label, x))

    return table[['start_time', 'length', 'label']]

def parse_interval_segments(day_segments):
    return day_segments

def parse_event_segments(day_segments):
    return day_segments

def parse_day_segments(day_segments_file, segments_type):
    # Add code to validate and parse frequencies, intervals, and events
    # Expected formats:
    # Frequency: label, length columns (e.g. my_prefix, 5) length has to be in minutes (int)
    # Interval: label, start, end columns (e.g. daily, 00:00, 23:59) start and end should be valid hours in 24 hour format
    # Event: label, timestamp, length, shift (e.g., survey1, 1532313215463, 60, -30), timestamp is a UNIX timestamp in ms (we could take a date time string instead), length is in minutes (int), shift is in minutes (+/-int) and is added/substracted from timestamp
    # Our output should have local_date, start_time, end_time, label. In the readable_datetime script, If local_date has the same value for all rows, every segment will be applied for all days, otherwise each segment will be applied only to its local_date
    day_segments = pd.read_csv(day_segments_file)

    if day_segments is None:
        message = 'The day segments file in [DAY_SEGMENTS][FILE] is None. Modify {}'.format(local_date)
        raise ValueError(message)

    if day_segments.shape[0] == 0:
        message = 'The day segments file in [DAY_SEGMENTS][FILE] is empty. Modify {}'.format(local_date)
        raise ValueError(message)

    if(segments_type not in ["FREQUENCY_EVERY_DAY", "INTERVAL_EVERY_DAY", "INTERVAL_FLEXIBLE_DAY"]):
        raise ValueError("[DAY_SEGMENTS][TYPE] can only be FREQUENCY_EVERY_DAY, INTERVAL_EVERY_DAY, or INTERVAL_FLEXIBLE_DAY")
    
    if(segments_type == "FREQUENCY_EVERY_DAY" and is_valid_frequency_segments(day_segments, day_segments_file)):
        day_segments = parse_frequency_segments(day_segments)
    elif(segments_type == "INTERVAL_EVERY_DAY" and is_valid_interval_segments(day_segments, day_segments_file)):
        day_segments = parse_interval_segments(day_segments)
    elif(segments_type == "INTERVAL_FLEXIBLE_DAY" and is_valid_event_segments(day_segments, day_segments_file)):
        day_segments = parse_event_segments(day_segments)
    else:
        raise ValueError("{} does not have a format compatible with frequency, interval or event day segments. Please refer to [LINK]".format(day_segments_file))
    return day_segments

final_day_segments = parse_day_segments(snakemake.input[0], snakemake.params["day_segments_type"])
final_day_segments.to_csv(snakemake.output["segments_file"], index=False)
pd.DataFrame({"label" : final_day_segments["label"].unique()}).to_csv(snakemake.output["segments_labels_file"], index=False)