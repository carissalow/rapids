import pandas as pd
import numpy as np

# def rapids_features(conversation_data, day_segment, requested_features,recordingMinutes,pausedMinutes,expectedMinutes):
def rapids_features(conversation_data, day_segment, provider, filter_data_by_segment, *args, **kwargs):

    requested_features = provider["FEATURES"]
    recordingMinutes = provider["RECORDING_MINUTES"]
    pausedMinutes = provider["PAUSED_MINUTES"]
    expectedMinutes =  1440 / (recordingMinutes + pausedMinutes)  

    # name of the features this function can compute
    base_features_names = ["minutessilence", "minutesnoise", "minutesvoice", "minutesunknown","sumconversationduration","avgconversationduration",
    "sdconversationduration","minconversationduration","maxconversationduration","timefirstconversation","timelastconversation","noisesumenergy",
    "noiseavgenergy","noisesdenergy","noiseminenergy","noisemaxenergy","voicesumenergy",
    "voiceavgenergy","voicesdenergy","voiceminenergy","voicemaxenergy","silencesensedfraction","noisesensedfraction",
    "voicesensedfraction","unknownsensedfraction","silenceexpectedfraction","noiseexpectedfraction","voiceexpectedfraction",
    "unknownexpectedfraction","countconversation"]

    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))

    conversation_features = pd.DataFrame(columns=["local_segment"] + ["conversation_rapids"  + "_" + x for x in features_to_compute])
    if not conversation_data.empty:
        conversation_data = filter_data_by_segment(conversation_data, day_segment)

        if not conversation_data.empty:
            conversation_features = pd.DataFrame()

            conversation_data = conversation_data.drop_duplicates(subset=["local_date", "local_time"], keep="first")

            if "minutessilence" in features_to_compute:
                conversation_features["conversation_rapids" + "_minutessilence"] = conversation_data[conversation_data['inference']==0].groupby(["local_segment"])['inference'].count()/60
            
            if "minutesnoise" in features_to_compute:
                conversation_features["conversation_rapids" + "_minutesnoise"] = conversation_data[conversation_data['inference']==1].groupby(["local_segment"])['inference'].count()/60

            if "minutesvoice" in features_to_compute:
                conversation_features["conversation_rapids" + "_minutesvoice"] = conversation_data[conversation_data['inference']==2].groupby(["local_segment"])['inference'].count()/60

            if "minutesunknown" in features_to_compute:
                conversation_features["conversation_rapids" + "_minutesunknown"] = conversation_data[conversation_data['inference']==3].groupby(["local_segment"])['inference'].count()/60

            if "countconversation" in features_to_compute:
                conversation_features["conversation_rapids" + "_countconversation"] = conversation_data[conversation_data["double_convo_start"] > 0].groupby(["local_segment"])['double_convo_start'].nunique()

            conv_duration = (conversation_data['double_convo_end']/1000 - conversation_data['double_convo_start']/1000)/60
            conversation_data = conversation_data.assign(conv_duration = conv_duration.values)
            
            conv_totalDuration = conversation_data[(conversation_data['inference'] >= 0) & (conversation_data['inference'] < 4)].groupby(["local_segment"])['inference'].count()/60 
            
            if "silencesensedfraction" in features_to_compute:
                conversation_features["conversation_rapids" + "_silencesensedfraction"] = (conversation_data[conversation_data['inference']==0].groupby(["local_segment"])['inference'].count()/60)/ conv_totalDuration

            if "noisesensedfraction" in features_to_compute:
                conversation_features["conversation_rapids" + "_noisesensedfraction"] = (conversation_data[conversation_data['inference']==1].groupby(["local_segment"])['inference'].count()/60)/ conv_totalDuration

            if "voicesensedfraction" in features_to_compute:
                conversation_features["conversation_rapids" + "_voicesensedfraction"] = (conversation_data[conversation_data['inference']==2].groupby(["local_segment"])['inference'].count()/60)/ conv_totalDuration

            if "unknownsensedfraction" in features_to_compute:
                conversation_features["conversation_rapids" + "_unknownsensedfraction"] = (conversation_data[conversation_data['inference']==3].groupby(["local_segment"])['inference'].count()/60)/ conv_totalDuration

            if "silenceexpectedfraction" in features_to_compute:
                conversation_features["conversation_rapids" + "_silenceexpectedfraction"] = (conversation_data[conversation_data['inference']==0].groupby(["local_segment"])['inference'].count()/60)/ expectedMinutes

            if "noiseexpectedfraction" in features_to_compute:
                conversation_features["conversation_rapids" + "_noiseexpectedfraction"] = (conversation_data[conversation_data['inference']==1].groupby(["local_segment"])['inference'].count()/60)/ expectedMinutes

            if "voiceexpectedfraction" in features_to_compute:
                conversation_features["conversation_rapids" + "_voiceexpectedfraction"] = (conversation_data[conversation_data['inference']==2].groupby(["local_segment"])['inference'].count()/60)/ expectedMinutes

            if "unknownexpectedfraction" in features_to_compute:
                conversation_features["conversation_rapids" + "_unknownexpectedfraction"] = (conversation_data[conversation_data['inference']==3].groupby(["local_segment"])['inference'].count()/60)/ expectedMinutes
                
            if "sumconversationduration" in features_to_compute:
                conversation_features["conversation_rapids" + "_sumconversationduration"] = conversation_data.groupby(["local_segment"])["conv_duration"].sum()

            if "avgconversationduration" in features_to_compute:
                conversation_features["conversation_rapids" + "_avgconversationduration"] = conversation_data[conversation_data["conv_duration"] > 0].groupby(["local_segment"])["conv_duration"].mean()

            if "sdconversationduration" in features_to_compute:
                conversation_features["conversation_rapids" + "_sdconversationduration"] = conversation_data[conversation_data["conv_duration"] > 0].groupby(["local_segment"])["conv_duration"].std()

            if "minconversationduration" in features_to_compute:
                conversation_features["conversation_rapids" + "_minconversationduration"] = conversation_data[conversation_data["conv_duration"] > 0].groupby(["local_segment"])["conv_duration"].min()

            if "maxconversationduration" in features_to_compute:
                conversation_features["conversation_rapids" + "_maxconversationduration"] = conversation_data.groupby(["local_segment"])["conv_duration"].max()

            if "timefirstconversation" in features_to_compute:
                timestampsLastConversation = conversation_data[conversation_data["double_convo_start"] > 0].groupby(["local_segment"])['timestamp'].min()
                if len(list(timestampsLastConversation.index)) > 0:
                    for date in list(timestampsLastConversation.index):
                        lastimestamp =  timestampsLastConversation.loc[date]
                        lasttime = (conversation_data.query('timestamp == @lastimestamp', inplace = False))['local_time'].iat[0]
                        conversation_features.loc[date,"conversation_rapids" + "_timefirstconversation"] = int(lasttime.split(':')[0])*60 + int(lasttime.split(':')[1])
                else:
                    conversation_features["conversation_rapids" + "_timefirstconversation"] = np.nan

            if "timelastconversation" in features_to_compute:
                timestampsLastConversation = conversation_data[conversation_data["double_convo_start"] > 0].groupby(["local_segment"])['timestamp'].max()
                if len(list(timestampsLastConversation.index)) > 0:
                    for date in list(timestampsLastConversation.index):
                        lastimestamp =  timestampsLastConversation.loc[date]
                        lasttime = (conversation_data.query('timestamp == @lastimestamp', inplace = False))['local_time'].iat[0]
                        conversation_features.loc[date,"conversation_rapids" + "_timelastconversation"] = int(lasttime.split(':')[0])*60 + int(lasttime.split(':')[1])
                else:
                    conversation_features["conversation_rapids" + "_timelastconversation"] = np.nan
            
            if "noisesumenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_noisesumenergy"] = conversation_data[conversation_data['inference']==1].groupby(["local_segment"])["double_energy"].sum()

            if "noiseavgenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_noiseavgenergy"] = conversation_data[conversation_data['inference']==1].groupby(["local_segment"])["double_energy"].mean()

            if "noisesdenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_noisesdenergy"] = conversation_data[conversation_data['inference']==1].groupby(["local_segment"])["double_energy"].std()

            if "noiseminenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_noiseminenergy"] = conversation_data[conversation_data['inference']==1].groupby(["local_segment"])["double_energy"].min()

            if "noisemaxenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_noisemaxenergy"] = conversation_data[conversation_data['inference']==1].groupby(["local_segment"])["double_energy"].max()

            if "voicesumenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_voicesumenergy"] = conversation_data[conversation_data['inference']==2].groupby(["local_segment"])["double_energy"].sum()

            if "voiceavgenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_voiceavgenergy"] = conversation_data[conversation_data['inference']==2].groupby(["local_segment"])["double_energy"].mean()

            if "voicesdenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_voicesdenergy"] = conversation_data[conversation_data['inference']==2].groupby(["local_segment"])["double_energy"].std()

            if "voiceminenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_voiceminenergy"] = conversation_data[conversation_data['inference']==2].groupby(["local_segment"])["double_energy"].min()

            if "voicemaxenergy" in features_to_compute:
                conversation_features["conversation_rapids" + "_voicemaxenergy"] = conversation_data[conversation_data['inference']==2].groupby(["local_segment"])["double_energy"].max()


            conversation_features = conversation_features.reset_index()

    return conversation_features