# Change Log

## Next release
- Add support for Empatica devices (all sensors)
- Add logo
- Move Citation page to the Setup section
- Add `config.yaml` validation schema and documentation.
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