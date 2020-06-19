import pandas as pd
from conversation.conversation_base import base_conversation_features

conversation_data = pd.read_csv(snakemake.input[0], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
requested_features = snakemake.params["features"]
recordingMinutes = snakemake.params["recordingMinutes"]
pausedMinutes = snakemake.params["pausedMinutes"]
expectedMinutes =  1440 / (recordingMinutes + pausedMinutes)  
conversation_features = pd.DataFrame(columns=["local_date"])

conversation_features = conversation_features.merge(base_conversation_features(conversation_data, day_segment, requested_features,recordingMinutes,pausedMinutes,expectedMinutes), on="local_date", how="outer")
assert len(requested_features) + 1 == conversation_features.shape[1], "The number of features in the output dataframe (=" + str(conversation_features.shape[1]) + ") does not match the expected value (=" + str(len(requested_features)) + " + 1). Verify your conversation feature extraction functions"

conversation_features.to_csv(snakemake.output[0], index=False)