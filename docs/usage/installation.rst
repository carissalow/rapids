Installation
===============

This instructions have been tested on MacOS Cataline and Ubuntu 16.04. If you find a problem, please report it.

Install MySQL

- ``brew install mysql``
- ``brew services start mysql``

Install R:

- ``brew install r``

Install miniconda:

- ``brew cask install miniconda`` (for Mac)
	
Clone our repo:

- ``git clone https://github.com/JulioV/moshi-aware.git``

Create a python virtual environment:

- ``conda env create -f environment.yml``
- ``conda init zsh/bash``
- Restart terminal if necessary
- ``conda activate moshi-env``

Install r packages and virtual environment:

- ``snakemake packrat_install``
- ``snakemake packrat_init``
- ``snakemake packrat_restore``

Configure the participants to analyze:

- Add a file per participant and name it with the label for that participant, no extension is necessary: ``/data/external/pxx``
- Add a line with the device_id for that participant as it appears on the db
- Add a line with the mobile platform (Android, or iOS)

Configure the db connection:

- Create an empty .env file in the root folder
- Add and complete the following lines
    | ``[MY_GROUP_NAME]``
    | ``user=MyUSER``
    | ``password=MyPassword``
    | ``host=MyIP``
    | ``port=3306``
