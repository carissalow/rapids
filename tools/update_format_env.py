import configparser
import yaml
from pathlib import Path

if not Path('./.env').exists():
    print("There is no .env file in RAPIDS root directory, nothing to migrate")
else:
    config = configparser.ConfigParser()
    config.read('./.env')

    credentials = dict()
    for section in config.sections():
        credentials[section] = dict()
        for attribute in config[section]:
            credentials[section][attribute] = config[section][attribute]

    with open('./credentials.yaml', 'w') as f:
        data = yaml.dump(credentials, f)

    print("Migration complete. Your credentials stored in .env should now be in credentials.yaml")