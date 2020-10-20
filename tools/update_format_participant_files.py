#!/usr/bin/python

from pathlib import Path
import yaml, os
import sys
p = Path(r'data/external/').glob('*')
files = [x for x in p if x.is_file() and x.suffix == "" and "." not in x.stem]
for file in files:
    reader = open(file, 'r') 
    phone = {"DEVICES_IDS" :"", "PLATFORMS" :"", "LABEL" :"", "START_DATE" :"", "END_DATE" :""}
    lines = reader.read().splitlines()
    if(len(lines) >=1 and len(lines[0]) > 0):
        phone["DEVICE_IDS"] = lines[0]
    if(len(lines) >=2 and len(lines[1]) > 0):
        phone["PLATFORMS"] = lines[1]
    if(len(lines) >=3 and len(lines[2]) > 0):
        phone["LABEL"] = lines[2]
    if(len(lines) >=4 and len(lines[3]) > 0):
        phone["START_DATE"] = lines[3].split(",")[0]
        phone["END_DATE"] = lines[3].split(",")[1]
    new_participant_file = Path(r'data/external/participant_files/') / (file.stem +  ".yaml")
    os.makedirs(os.path.dirname(new_participant_file), exist_ok=True)
    with open(new_participant_file, 'w') as writer:
        writer.write("PHONE:\n")
        writer.write("  DEVICE_IDS: [{}]\n".format(phone["DEVICE_IDS"]))
        writer.write("  PLATFORMS: [{}]\n".format(phone["PLATFORMS"]))
        writer.write("  LABEL: {}\n".format(phone["LABEL"]))
        writer.write("  START_DATE: {}\n".format(phone["START_DATE"]))
        writer.write("  END_DATE: {}\n".format(phone["END_DATE"]))

        writer.write("FITBIT:\n")
        writer.write("  DEVICE_IDS: [{}]\n".format(phone["DEVICE_IDS"]))
        writer.write("  LABEL: {}\n".format(phone["LABEL"]))
        writer.write("  START_DATE: {}\n".format(phone["START_DATE"]))
        writer.write("  END_DATE: {}\n".format(phone["END_DATE"]))
print("Processed files:")    
print(list(map(str, files)))