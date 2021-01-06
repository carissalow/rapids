# Change Log

## Next version v0.3.2
- Update docker and linux instructions to use RSPM binary repo for for faster installation
- Update CI to create a release on a tagged push that passes the tests
- Clarify in DB credential configuration that we only support MySQL
- Add Windows installation instructions
- Fix bugs in the create_participants_file script
- Fix bugs in Fitbit data parsing.
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