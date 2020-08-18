import pandas as pd
import configparser
import subprocess
import os

# read database credentials
group = snakemake.params["group"]
config = configparser.ConfigParser()
config.read(snakemake.input["db_credentials"])

# bash command to create table and restore tables from sql file
checkdb_cmd = "mysql -h " + config[group]["host"] + " -u " + config[group]["user"] + " -p" + config[group]["password"] + " -e use " + config[group]["database"]
create_cmd = "mysql -h " + config[group]["host"] + " -u " + config[group]["user"] + " -p" + config[group]["password"] + " -e \"CREATE DATABASE IF NOT EXISTS " + config[group]["database"] + ";\""
restore_cmd = "mysql -h " + config[group]["host"] + " -u " + config[group]["user"] + " -p" + config[group]["password"] + " " + config[group]["database"] + " < data/external/rapids_example.sql"

try:
    print("Checking if " + config[group]["database"] + " database exists")
    subprocess.run(checkdb_cmd.split(), check = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError:
    print(config[group]["database"] + " database does not exist")
    print("Creating " + config[group]["database"] + " database")
    os.system(create_cmd)
    print(config[group]["database"] + " database created")
    print("Restoring rapids_example.sql")
    os.system(restore_cmd)
    print("rapids_example.sql restored in " + config[group]["database"] + " database")
else:
    raise ValueError(config[group]["database"] + " DB already exists")
