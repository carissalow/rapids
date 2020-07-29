Manage virtual environments
=============================

**Add new packages**

Try to install any new package using `conda install my_package`. If a package is not available in one of conda's channels you can install it with pip but make sure your virtual environment is active.

**Update your conda environment.yaml**

After installing a new package you can use the following command in your terminal to update your ``environment.yaml`` before publishing your pipeline. Note that we ignore the package version for ``libfortran`` to keep compatibility with Linux:

    ``conda env export --no-builds | sed 's/^.*libgfortran.*$/  - libgfortran/' >  environment.yml``

**Update and prune your conda environment from a environment.yaml file**

Execute the following command in your terminal. See https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#updating-an-environment

    ``conda env update --prefix ./env --file environment.yml  --prune``