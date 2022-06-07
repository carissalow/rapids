import pandas as pd
import numpy as np

def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    keyboard_data = pd.read_csv(sensor_data_files["sensor_data"], dtype={"current_text":object, "before_text":object})
    typing_session_duration = provider["TYPING_SESSION_DURATION"] * 1000 # convert seconds to milliseconds
    requested_features = provider["FEATURES"]

    # name of the features this function can compute
    base_features_names = ["sessioncount","averageinterkeydelay","averagesessionlength","changeintextlengthlessthanminusone","changeintextlengthequaltominusone","changeintextlengthequaltoone","changeintextlengthmorethanone","maxtextlength","lastmessagelength","totalkeyboardtouches"]

    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))

    keyboard_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    if not keyboard_data.empty:
        keyboard_data = filter_data_by_segment(keyboard_data, time_segment)
        
        if not keyboard_data.empty:
            keyboard_features = pd.DataFrame()

            keyboard_data["keyboardStrokeDuration"] = keyboard_data['timestamp'].shift(-1) - keyboard_data['timestamp']
            keyboard_data["sessionStart"]           = keyboard_data["keyboardStrokeDuration"].apply(lambda x: 1 if x >= typing_session_duration else 0)
            keyboard_data["hasAppChanged"]          = keyboard_data["package_name"].shift(1, fill_value=keyboard_data["package_name"].head(1)) != keyboard_data["package_name"]
            keyboard_data["sessionStart"]           = keyboard_data[['hasAppChanged','sessionStart']].apply(lambda row: 1 if row['hasAppChanged'] == True else row['sessionStart'],axis=1)
            keyboard_data["sessionNumber"]          = keyboard_data["sessionStart"].cumsum()
            keyboard_data["changeInText"]           = keyboard_data["current_text"].str.len() - 2 - keyboard_data['before_text'].str.len()
            keyboard_data["currentTextLength"]      = keyboard_data["current_text"].str.len() - 2

            if "sessioncount" in features_to_compute:
                keyboard_features['sessioncount'] = keyboard_data.groupby(['local_segment'])['sessionStart'].sum()

            if "averagesessionlength" in features_to_compute:
                keyboard_features["averagesessionlength"] = keyboard_data[keyboard_data['sessionStart'] == 0].groupby(['local_segment','sessionNumber'])['keyboardStrokeDuration'].sum().reset_index().groupby(['local_segment'])['keyboardStrokeDuration'].mean()

            if "averageinterkeydelay" in features_to_compute:
                keyboard_features['averageinterkeydelay'] = keyboard_data[keyboard_data['sessionStart'] == 0].groupby(['local_segment','sessionNumber'])['keyboardStrokeDuration'].mean().reset_index().groupby(['local_segment'])['keyboardStrokeDuration'].mean()

            if "changeintextlengthlessthanminusone" in features_to_compute:
                keyboard_features['changeintextlengthlessthanminusone'] = keyboard_data[keyboard_data.changeInText < -1].groupby(['local_segment','sessionNumber'])['changeInText'].count().reset_index().groupby(['local_segment'])['changeInText'].count()

            if "changeintextlengthequaltominusone" in features_to_compute:
                keyboard_features['changeintextlengthequaltominusone'] = keyboard_data[keyboard_data.changeInText == -1].groupby(['local_segment','sessionNumber'])['changeInText'].count().reset_index().groupby(['local_segment'])['changeInText'].count()

            if "changeintextlengthequaltoone" in features_to_compute:
                keyboard_features['changeintextlengthequaltoone'] = keyboard_data[keyboard_data.changeInText == 1].groupby(['local_segment','sessionNumber'])['changeInText'].count().reset_index().groupby(['local_segment'])['changeInText'].count()

            if "changeintextlengthmorethanone" in features_to_compute:
                keyboard_features['changeintextlengthmorethanone'] = keyboard_data[keyboard_data.changeInText > 1].groupby(['local_segment','sessionNumber'])['changeInText'].count().reset_index().groupby(['local_segment'])['changeInText'].count()

            if "maxtextlength" in features_to_compute:
                keyboard_features["maxtextlength"] = keyboard_data[keyboard_data.currentTextLength > 0].groupby(['local_segment','sessionNumber'])['currentTextLength'].max().reset_index().groupby(['local_segment'])['currentTextLength'].mean()

            if "lastmessagelength" in features_to_compute:
                keyboard_data_copy = keyboard_data[['local_segment','sessionNumber','currentTextLength']].copy()
                keyboard_data_copy = keyboard_data_copy.drop_duplicates(subset = ["sessionNumber"],keep="last")
                keyboard_features["lastmessagelength"] = keyboard_data_copy[keyboard_data_copy.currentTextLength > 0].groupby(['local_segment','sessionNumber'])['currentTextLength'].mean().reset_index().groupby(['local_segment'])['currentTextLength'].mean()
            
            if "totalkeyboardtouches" in features_to_compute:
                keyboard_features["totalkeyboardtouches"] = keyboard_data.groupby(['local_segment','sessionNumber'])['is_password'].count().reset_index().groupby(['local_segment'])['is_password'].mean()
            
            keyboard_features.fillna(value={"sessioncount": 0, "averagesessionlength": 0, "changeintextlengthlessthanminusone": 0, "changeintextlengthequaltominusone": 0, "changeintextlengthequaltoone": 0, "changeintextlengthmorethanone": 0, "maxtextlength": 0, "totalkeyboardtouches": 0}, inplace=True)
            keyboard_features = keyboard_features.reset_index()

    return keyboard_features
