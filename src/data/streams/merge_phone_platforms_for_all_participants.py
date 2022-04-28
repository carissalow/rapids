import pandas as pd


platforms_of_all_participants = pd.DataFrame()
for platforms_file_path in snakemake.input["platforms_files"]:
    platforms_per_participant = pd.read_csv(platforms_file_path)
    platforms_per_participant.insert(0, "pid", platforms_file_path.split("/")[3])
    platforms_of_all_participants = pd.concat([platforms_of_all_participants, platforms_per_participant], axis=0)

platforms_of_all_participants.to_csv(snakemake.output[0], index=False)
