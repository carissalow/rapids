import pandas as pd

def is_valid_frequency_segments(day_segments):
    """
    returns true if day_segment has the expected structure for generating frequency segments;
    raises ValueError exception otherwise.
    """
    if day_segments is None:
        message = 'Table of frequency segmentation info is None. ' \
                  'Check the file under DAY_SEGMENTS in config.yaml'
        raise ValueError(message)

    if day_segments.shape[0] == 0:
        message = 'Table of frequency segmentation info is empty. ' \
                  'Check the file under DAY_SEGMENTS in config.yaml'
        raise ValueError(message)
    if day_segments.shape[0] > 1:
        message = 'Table of frequency segmentation info provides multiple specification but only one is allowed. ' \
                  'Check the file under DAY_SEGMENTS in config.yaml'
        raise ValueError(message)

    if 'length' not in day_segments.columns:
        message = 'Table of frequency segmentation info must provide segment length. ' \
                  'Check the file under DAY_SEGMENTS in config.yaml'
        raise ValueError(message)
    if 'label' not in day_segments.columns:
        message = 'Table of frequency segmentation info must provide segment label. ' \
                  'Check the file under DAY_SEGMENTS in config.yaml'
        raise ValueError(message)

    if not pd.api.types.is_integer_dtype(day_segments.dtypes['length']):
        message = 'Only integer segment length is allowed in the table of frequency segmentation; ' \
                  'found {}. Check the file under DAY_SEGMENTS in config.yaml'.format(day_segments.dtypes['length'])
        raise ValueError(message)

    if day_segments.iloc[0].loc['length'] < 0:
        message = 'Only positive integer segment length is allowed in the table of frequency segmentation; ' \
                  'found {}. Check the file under DAY_SEGMENTS in config.yaml'.format(day_segments.iloc[0].loc['length'])
        raise ValueError(message)
    if day_segments.iloc[0].loc['length'] >= 1440:
        message = 'Segment length in the table of frequency segmentation should be shorter than a day (in minutes); ' \
                  'found {}. Check the file under DAY_SEGMENTS in config.yaml'.format(day_segments.iloc[0].loc['length'])
        raise ValueError(message)

    return True

def is_valid_interval_segments(day_segments):
    return True

def is_valid_event_segments(day_segments):
    return False


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
    table['end_time'] = table['start_time'].shift(-1)
    table = table.iloc[:-1, :]

    label = day_segments.loc[0, 'label']
    table['label'] = range(0, table.shape[0])
    table['label'] = table['label'].apply(lambda x: '{}_{:04}'.format(label, x))

    table['local_date'] = None

    return table[['local_date', 'start_time', 'end_time', 'label']]

def parse_interval_segments(day_segments):
    day_segments["local_date"] = 1
    day_segments = day_segments.rename(columns={"start": "start_time", "end":"end_time"})
    return day_segments

def parse_event_segments(day_segments):
    return day_segments

def parse_day_segments(day_segments_file):
    # Add code to validate and parse frequencies, intervals, and events
    # Expected formats:
    # Frequency: label, length columns (e.g. my_prefix, 5) length has to be in minutes (int)
    # Interval: label, start, end columns (e.g. daily, 00:00, 23:59) start and end should be valid hours in 24 hour format
    # Event: label, timestamp, length, shift (e.g., survey1, 1532313215463, 60, -30), timestamp is a UNIX timestamp in ms (we could take a date time string instead), length is in minutes (int), shift is in minutes (+/-int) and is added/substracted from timestamp
    # Our output should have local_date, start_time, end_time, label. In the readable_datetime script, If local_date has the same value for all rows, every segment will be applied for all days, otherwise each segment will be applied only to its local_date
    day_segments = pd.read_csv(day_segments_file)

    if(is_valid_frequency_segments(day_segments)):
        day_segments = parse_frequency_segments(day_segments)
    elif(is_valid_interval_segments(day_segments)):
        day_segments = parse_interval_segments(day_segments)
    elif(is_valid_event_segments(day_segments)):
        day_segments = parse_event_segments(day_segments)
    else:
        raise ValueError("{} does not have a format compatible with frequency, interval or event day segments. Please refer to [LINK]".format(day_segments_file))
    return day_segments

day_segments = parse_day_segments(snakemake.input[0])
day_segments.to_csv(snakemake.output["segments_file"], index=False)