import pandas as pd

def unify_accelerometer_data(data):
    """
    Modify accelerometer data to convert readings from G-forces to m/s².

    Parameters:
    data (DataFrame): A pandas DataFrame containing accelerometer data.

    Returns:
    DataFrame: Modified accelerometer data with values in m/s².
    """
    GRAVITY_CONSTANT = 9.8
    data[['double_values_0', 'double_values_1', 'double_values_2']] *= GRAVITY_CONSTANT
    return data

def main(data, stream_parameters=None):
    """
    Main function to process the data using unify_accelerometer_data.

    Parameters:
    data (DataFrame): Input DataFrame containing accelerometer data.
    stream_parameters (dict): Additional parameters if needed (optional).

    Returns:
    DataFrame: Processed accelerometer data with values in m/s².
    """
    return unify_accelerometer_data(data)
