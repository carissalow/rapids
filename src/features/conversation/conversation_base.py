import pandas as pd

def base_conversation_features(conversation_data, day_segment, requested_features,recordingMinutes,pausedMinutes,expectedMinutes):
    # name of the features this function can compute
    base_features_names = ["minutessilence", "minutesnoise", "minutesvoice", "minutesunknown","sumconversationduration","avgconversationduration",
    "sdconversationduration","minconversationduration","maxconversationduration","timefirstconversation","timelastconversation","sumenergy",
    "avgenergy","sdenergy","minenergy","maxenergy","silencesensedfraction","noisesensedfraction",
    "voicesensedfraction","unknownsensedfraction","silenceexpectedfraction","noiseexpectedfraction","voiceexpectedfraction",
    "unknownexpectedfraction"]

    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))


    if conversation_data.empty:
        conversation_features = pd.DataFrame(columns=["local_date"] + ["conversation_" + day_segment + "_" + x for x in features_to_compute])
    else:
        if day_segment != "daily":
            conversation_data = conversation_data[conversation_data["local_day_segment"] == day_segment]

        if conversation_data.empty:
            conversation_features = pd.DataFrame(columns=["local_date"] + ["conversation_" + day_segment + "_" + x for x in features_to_compute])
        else:
            conversation_features = pd.DataFrame()

            if "minutessilence" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_minutessilence"] = conversation_data[conversation_data["inference"]==0].groupby(["local_date"])["inference"].count()
            
            if "minutesnoise" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_minutesnoise"] = conversation_data[conversation_data["inference"]==1].groupby(["local_date"])["inference"].count()

            if "minutesvoice" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_minutesvoice"] = conversation_data[conversation_data["inference"]==2].groupby(["local_date"])["inference"].count()

            if "minutesunknown" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_minutesunknown"] = conversation_data[conversation_data["inference"]==3].groupby(["local_date"])["inference"].count()

            conv_duration = conversation_data["double_convo_end"] - conversation_data["double_convo_start"]
            conversation_data = conversation_data.assign(conv_duration = conv_duration.values)

            conversation_data["total_duration"] = conversation_data[conversation_data["inference"]==0].groupby(["local_date"])["inference"].count() + conversation_data[conversation_data["inference"]==1].groupby(["local_date"])["inference"].count() + conversation_data[conversation_data["inference"]==2].groupby(["local_date"])["inference"].count() + conversation_data[conversation_data["inference"]==3].groupby(["local_date"])["inference"].count()
            
            if "silencesensedfraction" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_silencesensedfraction"] = conversation_data[conversation_data["inference"]==0].groupby(["local_date"])["inference"].count()/ conversation_data["total_duration"]

            if "noisesensedfraction" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_noisesensedfraction"] = conversation_data[conversation_data["inference"]==1].groupby(["local_date"])["inference"].count()/ conversation_data["total_duration"]

            if "voicesensedfraction" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_voicesensedfraction"] = conversation_data[conversation_data["inference"]==2].groupby(["local_date"])["inference"].count()/ conversation_data["total_duration"]

            if "unknownsensedfraction" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_unknownsensedfraction"] = conversation_data[conversation_data["inference"]==3].groupby(["local_date"])["inference"].count()/ conversation_data["total_duration"]    

            if "silenceexpectedfraction" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_silenceexpectedfraction"] = conversation_data[conversation_data["inference"]==0].groupby(["local_date"])["inference"].count()/ expectedMinutes

            if "noiseexpectedfraction" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_noiseexpectedfraction"] = conversation_data[conversation_data["inference"]==1].groupby(["local_date"])["inference"].count()/ expectedMinutes

            if "voiceexpectedfraction" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_voiceexpectedfraction"] = conversation_data[conversation_data["inference"]==2].groupby(["local_date"])["inference"].count()/ expectedMinutes

            if "unknownexpectedfraction" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_unknownexpectedfraction"] = conversation_data[conversation_data["inference"]==3].groupby(["local_date"])["inference"].count()/ expectedMinutes
                
            if "sumconversationduration" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_sumconversationduration"] = conversation_data.groupby(["local_date"])["conv_duration"].sum()

            if "avgconversationduration" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_avgconversationduration"] = conversation_data.groupby(["local_date"])["conv_duration"].mean()

            if "sdconversationduration" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_sdconversationduration"] = conversation_data.groupby(["local_date"])["conv_duration"].std()

            if "minconversationduration" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_minconversationduration"] = conversation_data.groupby(["local_date"])["conv_duration"].min()

            if "maxconversationduration" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_maxconversationduration"] = conversation_data.groupby(["local_date"])["conv_duration"].max()

            if "timefirstconversation" in features_to_compute:
                timeFirstConversation = conversation_data[conversation_data["double_convo_start"] > 0].groupby(["local_date"])[["double_convo_start","local_hour","local_minute"]].min()
                if "local_hour" in timeFirstConversation.columns:
                    conversation_features["conversation_" + day_segment + "_timefirstconversation"] = timeFirstConversation["local_hour"]*60 + timeFirstConversation["local_minute"]
                else:
                    conversation_features["conversation_" + day_segment + "_timefirstconversation"] = 0

            if "timelastconversation" in features_to_compute:
                timeLastConversation = conversation_data[conversation_data["double_convo_start"] > 0].groupby(["local_date"])[["double_convo_start","local_hour","local_minute"]].max()
                if "local_hour" in timeLastConversation:
                    conversation_features["conversation_" + day_segment + "_timelastconversation"] = timeLastConversation["local_hour"]*60 + timeLastConversation["local_minute"]
                else:
                    conversation_features["conversation_" + day_segment + "_timelastconversation"] = 0

            if "sumenergy" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_sumenergy"] = conversation_data.groupby(["local_date"])["double_energy"].sum()

            if "avgenergy" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_avgenergy"] = conversation_data.groupby(["local_date"])["double_energy"].mean()

            if "sdenergy" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_sdenergy"] = conversation_data.groupby(["local_date"])["double_energy"].std()

            if "minenergy" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_minenergy"] = conversation_data.groupby(["local_date"])["double_energy"].min()

            if "maxenergy" in features_to_compute:
                conversation_features["conversation_" + day_segment + "_maxenergy"] = conversation_data.groupby(["local_date"])["double_energy"].max()


            conversation_features = conversation_features.reset_index()

    return conversation_features