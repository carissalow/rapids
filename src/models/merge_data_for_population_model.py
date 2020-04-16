import pandas as pd

data_all_participants = pd.DataFrame()
for data_file in snakemake.input["data_files"]:
    data_single_participant = pd.read_csv(data_file)
    data_all_participants = pd.concat([data_all_participants, data_single_participant], axis=0)

data_all_participants.to_csv(snakemake.output[0], index=False)
