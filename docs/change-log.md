# Change Log
## v1.2.1 (in progress)
- Fix PID matching when joining data from multiple participants. We can handle PIDS with an arbitrary format.
- Fix bug that did not correctly parse participants with more than 2 phones or more than 1 wearable
- New keyboard features
- Add the `EXCLUDE_SLEEP` module for steps intraday features
- Fix bug when no phone data yield is needed to process location data
- Remove location rows with the same timestamp based on their accuracy
## v1.2.0
- Sleep summary and intraday features are more consistent.
- Add wake and bedtime features for sleep summary data.
- Fix bugs with sleep PRICE features.
- Update home page
- Add contributing guide
## v1.1.1
- Fix length of periodic segments on days with DLS
- Fix crash when scraping data for an app that does not exist
- Add tests for phone screen data
## v1.1.0
- Add Fitbit calories intraday features
## v1.0.1
- Fix crash in `chunk_episodes` of `utils.py` for multi time zone data
- Fix crash in BT Doryab provider when the number of clusters is 2
- Fix Fitbit multi time zone inference from phone data (simplify)
- Fix missing columns when the input for phone data yield is empty
- Fix wrong date time labels for event segments for multi time zone data (all labels are computed based on a single tz)
- Fix periodic segment crash when there are no segments to assign (only affects wday, mday, qday, or yday) 
- Fix crash in Analysis Workflow with new suffix in segments' labels
## v1.0.0
- Add a new [Overview](../setup/overview/) page.
- You can [extend](../datastreams/add-new-data-streams/) RAPIDS with your own [data streams](../datastreams/data-streams-introduction/). Data streams are data collected with other sensing apps besides AWARE (like Beiwe, mindLAMP), and stored in other data containers (databases, files) besides MySQL.
- Support to analyze Empatica wearable data (thanks to Joe Kim and  Brinnae Bent from the [DBDP](https://dbdp.org/))
- Support to analyze AWARE data stored in [CSV files](../datastreams/aware-csv/) and [InfluxDB](../datastreams/aware-influxdb/) databases
- Support to analyze data collected over [multiple time zones](../setup/configuration/#multiple-timezones)
- Support for [sleep intraday features](../features/fitbit-sleep-intraday/) from the core team and also from the community (thanks to Stephen Price)
- Users can comment on the documentation (powered by utterances).
- `SCR_SCRIPT` and `SRC_LANGUAGE` are replaced by `SRC_SCRIPT`.
- Add RAPIDS new logo
- Move Citation and Minimal Example page to the Setup section
- Add `config.yaml` validation schema and documentation. Now it's more difficult to modify the `config.yaml` file with invalid values.
- Add new `time at home` Doryab location feature
- Add and home coordinates to the location data file so location providers can build features based on it.
- If you are migrating from RAPIDS 0.4.3 or older, check this [guide](../migrating-from-old-versions/#migrating-from-rapids-04x-or-older)
## v0.4.3
- Fix bug when any of the rows from any sensor do not belong a time segment
## v0.4.2
- Update battery testing
- Fix location processing bug when certain columns don't exist
- Fix HR intraday bug when minutesonZONE features were 0 
- Update FAQs
- Fix HR summary bug when restinghr=0 (ignore those rows)
- Fix ROG, location entropy and normalized entropy in Doryab location provider
- Remove sampling frequency dependance in Doryab location provider
- Update documentation of Doryab location provider
- Add new `FITBIT_DATA_YIELD` `RAPIDS` provider
- Deprecate Doryab circadian movement feature until it is fixed
## v0.4.1
- Fix bug when no error message was displayed for an empty `[PHONE_DATA_YIELD][SENSORS]` when resampling location data
## v0.4.0
- Add four new phone sensors that can be used for PHONE_DATA_YIELD
- Add code so new feature providers can be added for the new four sensors
- Add new clustering algorithm (OPTICS) for Doryab features
- Update default EPS parameter for Doryab location clustering
- Add clearer error message for invalid phone data yield sensors
- Add ALL_RESAMPLED flag and accuracy limit for location features
- Add FAQ about null characters in phone tables
- Reactivate light and wifi tests and update testing docs
- Fix bug when parsing Fitbit steps data
- Fix bugs when merging features from empty time segments
- Fix minor issues in the documentation
## v0.3.2
- Update docker and linux instructions to use RSPM binary repo for for faster installation
- Update CI to create a release on a tagged push that passes the tests
- Clarify in DB credential configuration that we only support MySQL
- Add Windows installation instructions
- Fix bugs in the create_participants_file script
- Fix bugs in Fitbit data parsing.
- Fixed Doryab location features context of clustering.
- Fixed the wrong shifting while calculating distance in Doryab location features.
- Refactored the haversine function
## v0.3.1
- Update installation docs for RAPIDS' docker container
- Fix example analysis use of accelerometer data in a plot
- Update FAQ
- Update minimal example documentation
- Minor doc updates
## v0.3.0
- Update R and Python virtual environments
- Add GH actions CI support for tests and docker
- Add release and test badges to README
## v0.2.6
- Fix old versions banner on nested pages
## v0.2.5
- Fix docs deploy typo
## v0.2.4
- Fix broken links in landing page and docs deploy
## v0.2.3
- Fix participant IDS in the example analysis workflow
## v0.2.2
- Fix readme link to docs
## v0.2.1
- FIx link to the most recent version in the old version banner

## v0.2.0
- Add new `PHONE_BLUETOOTH` `DORYAB` provider
- Deprecate `PHONE_BLUETOOTH` `RAPIDS` provider
- Fix bug in `filter_data_by_segment` for Python when dataset was empty
- Minor doc updates
- New FAQ item

## v0.1.0
- New and more consistent docs (this website). The [previous docs](https://rapidspitt.readthedocs.io/en/latest/) are marked as beta 
- Consolidate [configuration](../setup/configuration) instructions
- Flexible [time segments](../setup/configuration#time-segments)
- Simplify Fitbit behavioral feature extraction and [documentation](../features/fitbit-heartrate-summary)
- Sensor's configuration and output is more consistent
- Update [visualizations](../visualizations/data-quality-visualizations) to handle flexible day segments
- Create a RAPIDS [execution](../setup/execution) script that allows re-computation of the pipeline after configuration changes
- Add [citation](../citation) guide
- Update [virtual environment](../developers/virtual-environments) guide
- Update analysis workflow [example](../workflow-examples/analysis)
- Add a [Code of Conduct](../code_of_conduct)
- Update [Team](../team) page