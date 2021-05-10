# Contributing

Thank you for taking the time to contribute! 

All changes, small or big, are welcome, and regardless of who you are, we are always happy to work together to make your contribution as strong as possible. We follow the [Covenant Code of Conduct](../code_of_conduct), so we ask you to uphold it. Be kind to everyone in the community, and please report unacceptable behavior to moshiresearch@gmail.com.

## Questions, Feature Requests, and Discussions

Post any questions, feature requests, or discussions in our [GitHub Discussions tab](https://github.com/carissalow/rapids/discussions).

## Bug Reports

Report any bugs in our [GithHub issue tracker](https://github.com/carissalow/rapids/issues) keeping in mind to:

- Debug and simplify the problem to create a minimal example. For example, reduce the problem to a single participant, sensor, and a few rows of data.
- Provide a clear and succinct description of the problem (expected behavior vs. actual behavior).
- Attach your `config.yaml`, time segments file, and time zones file if appropriate.
- Attach test data if possible and any screenshots or extra resources that will help us debug the problem.
- Share the commit you are running: `git rev-parse --short HEAD`
- Share your OS version (e.g., Windows 10)
- Share the device/sensor you are processing (e.g., phone accelerometer)

## Documentation Contributions

If you want to fix a typo or any other minor changes, you can edit the file online by clicking on the pencil icon at the top right of any page and opening a pull request using [Github's website](https://docs.github.com/en/github/managing-files-in-a-repository/editing-files-in-your-repository)

If your changes are more complex, clone RAPIDS' repository, setup the dev environment for our documentation with this [tutorial](../developers/documentation), and submit any changes on a new *feature branch* following our [git flow](../developers/git-flow).

## Code Contributions

!!! hint "Hints for any code changes"
    - To submit any new code, use a new *feature branch* following our [git flow](../developers/git-flow).
    - If you neeed a new Python or R package in RAPIDS' virtual environments, follow this [tutorial](../developers/virtual-environments/)
    - If you need to change the `config.yaml` you will need to update its validation schema with this [tutorial](../developers/validation-schema-config/)

### New Data Streams

*New data containers.* If you want to process data from a device RAPIDS supports ([see this table](../datastreams/data-streams-introduction/)) but it's stored in a database engine or file type we don't support yet, [implement a new data stream container and format](../datastreams/add-new-data-streams/). You can copy and paste the `format.yaml` of one of the other streams of the device you are targeting.

*New sensing apps.* If you want to add support for new smartphone sensing apps like Beiwe, [implement a new data stream container and format](../datastreams/add-new-data-streams/).

*New wearable devices.* If you want to add support for a new wearable, open a [Github discussion](https://github.com/carissalow/rapids/discussions), so we can add the necessary initial configuration files and code.

### New Behavioral Features

If you want to add new [behavioral features](../features/feature-introduction/) for mobile sensors RAPIDS already supports, follow this [tutorial](../features/add-new-features/). A sensor is supported if it has a configuration section in `config.yaml`.

If you want to add new [behavioral features](../features/feature-introduction/) for mobile sensors RAPIDS does not support yet, open a [Github discussion](https://github.com/carissalow/rapids/discussions), so we can add the necessary initial configuration files and code.

### New Tests

If you want to add new tests for existent behavioral features, follow this [tutorial](../developers/testing).

### New Visualizations

Open a [Github discussion](https://github.com/carissalow/rapids/discussions), so we can add the necessary initial configuration files and code.