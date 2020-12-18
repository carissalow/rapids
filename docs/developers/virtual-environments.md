## Python Virtual Environment

### Add new packages
Try to install any new package using `conda install -c CHANNEL PACKAGE_NAME` (you can use `pip` if the package is only available there). Make sure your Python virtual environment is active (`conda activate YOUR_ENV`).

### Remove packages
Uninstall packages using the same manager you used to install them `conda remove PACKAGE_NAME` or `pip uninstall PACKAGE_NAME`

### Updating all packages
Make sure your Python virtual environment is active (`conda activate YOUR_ENV`), then run
```bash
conda update --all
```

### Update your conda `environment.yaml`
After installing or removing a package you can use the following command in your terminal to update your `environment.yaml` before publishing your pipeline. Note that we ignore the package version for `libfortran` and `mkl` to keep compatibility with Linux:
```bash
conda env export --no-builds | sed 's/^.*libgfortran.*$/  - libgfortran/' | sed 's/^.*mkl=.*$/  - mkl/' >  environment.yml
```

## R Virtual Environment

### Add new packages
1. Open your terminal and navigate to RAPIDS' root folder
2. Run `R` to open an R interactive session
3. Run `renv::install("PACKAGE_NAME")`

### Remove packages
1. Open your terminal and navigate to RAPIDS' root folder
2. Run `R` to open an R interactive session
3. Run `renv::remove("PACKAGE_NAME")`

### Updating all packages
1. Open your terminal and navigate to RAPIDS' root folder
2. Run `R` to open an R interactive session
3. Run `renv::update()`
### Update your R `renv.lock`
After installing or removing a package you can use the following command in your terminal to update your `renv.lock` before publishing your pipeline.

1. Open your terminal and navigate to RAPIDS' root folder
2. Run `R` to open an R interactive session
3. Run `renv::snapshot()` (renv will ask you to confirm any updates to this file)

