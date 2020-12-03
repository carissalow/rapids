# Migrating from RAPIDS beta

If you were relying on the [old docs](https://rapidspitt.readthedocs.io/en/latest/) and the most recent version of RAPIDS you downloaded is from or before [Oct 13, 2020](https://github.com/carissalow/rapids/commit/640890c7b49492d150accff5c87b1eb25bd97a49) you are using the beta version of RAPIDS.

You can start using the new RAPIDS (we are starting with `v0.1.0`) right away, just take into account the following:

1. [Install](setup/installation.md) a new copy of RAPIDS (the R and Python virtual environments didn't change so the cached versions will be reused)
      1. Make sure you don't skip a new Installation step to give execution permissions to the RAPIDS script: `chmod +x rapids`
2. Follow the new [Configuration](setup/configuration.md) guide.
      1. You can copy and paste your old `.env` file
      2. You can migrate your old participant files: 
      ```
      python tools/update_format_participant_files.py
      ```
3. You can proceed to reconfigure your `config.yaml`, its structure is more consistent and should be familiar to you.

!!! info
    If you have any questions reach out to us on [Slack](http://awareframework.com:3000/).