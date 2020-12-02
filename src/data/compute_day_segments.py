import pandas as pd
import warnings
import yaml

def is_valid_frequency_segments(day_segments, day_segments_file):
    """
    returns true if day_segment has the expected structure for generating frequency segments;
    raises ValueError exception otherwise.
    """
    
    valid_columns = ["label", "length"]
    if set(day_segments.columns) != set(valid_columns):
        error_message = 'The FREQUENCY day segments file in [DAY_SEGMENTS][FILE] must have two columns: label, and length ' \
                  'but instead we found {}. Modify {}'.format(list(day_segments.columns), day_segments_file)
        raise ValueError(error_message)

    if day_segments.shape[0] > 1:
        message = 'The FREQUENCY day segments file in [DAY_SEGMENTS][FILE] can only have 1 row.' \
                  'Modify {}'.format(day_segments_file)
        raise ValueError(message)

    if not pd.api.types.is_integer_dtype(day_segments.dtypes['length']):
        message = 'The column length in the FREQUENCY day segments file in [DAY_SEGMENTS][FILE] must be integer but instead is ' \
                  '{}. . This usually means that not all values in this column are formed by digits. Modify {}'.format(day_segments.dtypes['length'], day_segments_file)
        raise ValueError(message)

    if day_segments.iloc[0].loc['length'] < 0:
        message = 'The value in column length in the FREQUENCY day segments file in [DAY_SEGMENTS][FILE] must be positive but instead is  ' \
                  '{}. Modify {}'.format(day_segments.iloc[0].loc['length'], day_segments_file)
        raise ValueError(message)
    if day_segments.iloc[0].loc['length'] >= 1440:
        message = 'The column length in the FREQUENCY day segments file in [DAY_SEGMENTS][FILE] must be shorter than a day in minutes (1440) but instead is ' \
                  '{}. Modify {}'.format(day_segments.iloc[0].loc['length'], day_segments_file)
        raise ValueError(message)

    return True

def is_valid_periodic_segments(day_segments, day_segments_file):
    day_segments = day_segments.copy(deep=True)

    valid_columns = ["label", "start_time", "length", "repeats_on", "repeats_value"]
    if set(day_segments.columns) != set(valid_columns):
        error_message = 'The PERIODIC day segments file in [DAY_SEGMENTS][FILE] must have five columns: label, start_time, length, repeats_on, repeats_value ' \
                  'but instead we found {}. Modify {}'.format(list(day_segments.columns), day_segments_file)
        raise ValueError(error_message)

    valid_repeats_on = ["every_day", "wday", "mday", "qday", "yday"]
    if len(list(set(day_segments["repeats_on"]) - set(valid_repeats_on))) > 0:
        error_message = 'The column repeats_on in the PERIODIC day segments file in [DAY_SEGMENTS][FILE] can only accept: "every_day", "wday", "mday", "qday", or "yday" ' \
                  'but instead we found {}. Modify {}'.format(list(set(day_segments["repeats_on"])), day_segments_file)
        raise ValueError(error_message)

    if not pd.api.types.is_integer_dtype(day_segments.dtypes['repeats_value']):
        message = 'The column repeats_value in the PERIODIC day segments file in [DAY_SEGMENTS][FILE] must be integer but instead is ' \
                  '{}. . This usually means that not all values in this column are formed by digits. Modify {}'.format(day_segments.dtypes['repeats_value'], day_segments_file)
        raise ValueError(message)

    invalid_day_segments = day_segments.query("repeats_on == 'every_day' and repeats_value != 0")
    if invalid_day_segments.shape[0] > 0:
        message = 'Every row with repeats_on=every_day must have a repeats_value=0 in the PERIODIC day segments file in [DAY_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_day_segments["label"].to_numpy(), day_segments_file)
        raise ValueError(message)

    invalid_day_segments = day_segments.query("repeats_on == 'wday' and (repeats_value < 1 | repeats_value > 7)")
    if invalid_day_segments.shape[0] > 0:
        message = 'Every row with repeats_on=wday must have a repeats_value=[1,7] in the PERIODIC day segments file in [DAY_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_day_segments["label"].to_numpy(), day_segments_file)
        raise ValueError(message)

    invalid_day_segments = day_segments.query("repeats_on == 'mday' and (repeats_value < 1 | repeats_value > 31)")
    if invalid_day_segments.shape[0] > 0:
        message = 'Every row with repeats_on=mday must have a repeats_value=[1,31] in the PERIODIC day segments file in [DAY_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_day_segments["label"].to_numpy(), day_segments_file)
        raise ValueError(message)

    invalid_day_segments = day_segments.query("repeats_on == 'qday' and (repeats_value < 1 | repeats_value > 92)")
    if invalid_day_segments.shape[0] > 0:
        message = 'Every row with repeats_on=qday must have a repeats_value=[1,92] in the PERIODIC day segments file in [DAY_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_day_segments["label"].to_numpy(), day_segments_file)
        raise ValueError(message)

    invalid_day_segments = day_segments.query("repeats_on == 'yday' and (repeats_value < 1 | repeats_value > 366)")
    if invalid_day_segments.shape[0] > 0:
        message = 'Every row with repeats_on=yday must have a repeats_value=[1,366] in the PERIODIC day segments file in [DAY_SEGMENTS][FILE].' \
                  ' Modify row(s) of segment(s) {} of {}'.format(invalid_day_segments["label"].to_numpy(), day_segments_file)
        raise ValueError(message)

    try:
        day_segments["start_time"] = pd.to_datetime(day_segments["start_time"])
    except ValueError as err:
        raise ValueError("At least one start_time in the PERIODIC day segments file in [DAY_SEGMENTS][FILE] has an invalid format, it should be HH:MM:SS in 24hr clock({}). Modify {}".format(err, day_segments_file))

    if(day_segments.shape[0] != day_segments.drop_duplicates().shape[0]):
        error_message = 'The PERIODIC day segments file in [DAY_SEGMENTS][FILE] has two or more rows that are identical. ' \
                  'Modify {}'.format(day_segments_file)
        raise ValueError(error_message)
    
    duplicated_labels = day_segments[day_segments["label"].duplicated()]
    if(duplicated_labels.shape[0] > 0):
        error_message = 'Segements labels must be unique. The PERIODIC day segments file in [DAY_SEGMENTS][FILE] has {} row(s) with the same label {}. ' \
                  'Modify {}'.format(duplicated_labels.shape[0], duplicated_labels["label"].to_numpy(), day_segments_file)
        raise ValueError(error_message)
    
    # TODO Validate string format for lubridate
    
    return True

def is_valid_event_segments(day_segments, day_segments_file):
    day_segments = day_segments.copy(deep=True)

    valid_columns = ["label", "event_timestamp", "length", "shift", "shift_direction", "device_id"]
    if set(day_segments.columns) != set(valid_columns):
        error_message = 'The EVENT day segments file in [DAY_SEGMENTS][FILE] must have six columns: label, event_timestamp, length, shift, shift_direction and device_id ' \
                  'but instead we found {}. Modify {}'.format(list(day_segments.columns), day_segments_file)
        raise ValueError(error_message)

    if not pd.api.types.is_integer_dtype(day_segments.dtypes['event_timestamp']):
        message = 'The column event_timestamp in the EVENT day segments file in [DAY_SEGMENTS][FILE] must be integer but instead is ' \
                  '{}. This usually means that not all values in this column are formed by digits. Modify {}'.format(day_segments.dtypes['event_timestamp'], day_segments_file)
        raise ValueError(message)

    valid_shift_direction_values = [1, -1, 0]
    provided_values = day_segments["shift_direction"].unique()
    if len(list(set(provided_values) - set(valid_shift_direction_values))) > 0:
        error_message = 'The values of shift_direction column in the EVENT day segments file in [DAY_SEGMENTS][FILE] can only be 1, -1 or 0 ' \
                  'but instead we found {}. Modify {}'.format(provided_values, day_segments_file)
        raise ValueError(error_message)

    if(day_segments.shape[0] != day_segments.drop_duplicates().shape[0]):
        error_message = 'The EVENT day segments file in [DAY_SEGMENTS][FILE] has two or more rows that are identical. ' \
                  'Modify {}'.format(day_segments_file)
        raise ValueError(error_message)

    # TODO Validate string format for lubridate of length and shift
    # TODO validate unique labels per participant
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

def parse_periodic_segments(day_segments):
    day_segments.loc[day_segments["repeats_on"] == "every_day", "repeats_value"] = 0
    return day_segments

def parse_event_segments(day_segments, device_ids):
    return day_segments.query("device_id == @device_ids")

def parse_day_segments(day_segments_file, segments_type, device_ids):
    # Add code to validate and parse frequencies, intervals, and events
    # Expected formats:
    # Frequency: label, length columns (e.g. my_prefix, 5) length has to be in minutes (int)
    # Interval: label, start, end columns (e.g. daily, 00:00, 23:59) start and end should be valid hours in 24 hour format
    # Event: label, timestamp, length, shift (e.g., survey1, 1532313215463, 60, -30), timestamp is a UNIX timestamp in ms (we could take a date time string instead), length is in minutes (int), shift is in minutes (+/-int) and is added/substracted from timestamp
    # Our output should have local_date, start_time, end_time, label. In the readable_datetime script, If local_date has the same value for all rows, every segment will be applied for all days, otherwise each segment will be applied only to its local_date
    day_segments = pd.read_csv(day_segments_file)

    if day_segments is None:
        message = 'The day segments file in [DAY_SEGMENTS][FILE] is None. Modify {}'.format(day_segments_file)
        raise ValueError(message)

    if day_segments.shape[0] == 0:
        message = 'The day segments file in [DAY_SEGMENTS][FILE] is empty. Modify {}'.format(day_segments_file)
        raise ValueError(message)

    if(segments_type not in ["FREQUENCY", "PERIODIC", "EVENT"]):
        raise ValueError("[DAY_SEGMENTS][TYPE] can only be FREQUENCY, PERIODIC, or EVENT")
    
    if(segments_type == "FREQUENCY" and is_valid_frequency_segments(day_segments, day_segments_file)):
        day_segments = parse_frequency_segments(day_segments)
    elif(segments_type == "PERIODIC" and is_valid_periodic_segments(day_segments, day_segments_file)):
        day_segments = parse_periodic_segments(day_segments)
    elif(segments_type == "EVENT" and is_valid_event_segments(day_segments, day_segments_file)):
        day_segments = parse_event_segments(day_segments, device_ids)
    else:
        raise ValueError("{} does not have a format compatible with frequency, periodic or event day segments. Please refer to [LINK]".format(day_segments_file))
    return day_segments

participant_file = yaml.load(open(snakemake.input[1], 'r'), Loader=yaml.FullLoader)
device_ids = []
for key in participant_file.keys():
    if "DEVICE_IDS" in participant_file[key]:
        device_ids = device_ids + participant_file[key]["DEVICE_IDS"]

final_day_segments = parse_day_segments(snakemake.input[0], snakemake.params["day_segments_type"], device_ids)

if snakemake.params["day_segments_type"] == "EVENT" and final_day_segments.shape[0] == 0:
    warnings.warn("There are no event day segments for {}. Check your day segment file {}".format(snakemake.params["pid"], snakemake.input[0]))

final_day_segments.to_csv(snakemake.output["segments_file"], index=False)
pd.DataFrame({"label" : final_day_segments["label"].unique()}).to_csv(snakemake.output["segments_labels_file"], index=False)