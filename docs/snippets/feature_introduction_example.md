1. **Sensor section**

    Each sensor (accelerometer, screen, etc.) of every supported device (smartphone, Fitbit, etc.) has a section in the `config.yaml` with `parameters` and feature `PROVIDERS`.

2. **Sensor Parameters.** 
    
    Each sensor section has one or more parameters. These are parameters that affect different aspects of how the raw data is pulled, and processed.
    
    The `CONTAINER` parameter exists for every sensor, but some sensors will have extra parameters like [`[PHONE_LOCATIONS]`](../phone-locations/).
    
    We explain these parameters in a table at the top of each sensor documentation page.

3. **Sensor Providers**

    Each object in this list represents a feature `PROVIDER`. Each sensor can have zero, one, or more providers.
    
    A `PROVIDER` is a script that creates behavioral features for a specific sensor. Providers are created by the core RAPIDS team or by the community, which are named after its first author like [[PHONE_LOCATIONS][DORYAB]](../../features/phone-locations/#doryab-provider).

    In this example, there are two accelerometer feature providers `RAPIDS` and `PANDA`.

4. **`PROVIDER` Parameters**
    
    Each `PROVIDER` has parameters that affect the computation of the behavioral features it offers.
    
    These parameters include at least a `[COMPUTE]` flag that you switch to `True` to extract a provider's behavioral features. 

    We explain every provider's parameter in a table under the `Parameters description` heading on each provider documentation page.

5. **`PROVIDER` Features**

    Each `PROVIDER` offers a set of behavioral features.
    
    These features are grouped in an array for some providers, like those for `RAPIDS` provider. For others, they are grouped in a collection of arrays, like those for `PANDAS` provider.
    
    In either case, you can delete the features you are not interested in, and they will not be included in the sensor's output feature file. 

    We explain each behavioral feature in a table under the `Features description` heading on each provider documentation page.

6. **`PROVIDER` script**

    Each `PROVIDER` has a `SRC_SCRIPT` that points to the script implementing its behavioral features.

    It has to be a relative path from RAPIDS' root folder and the script's parent folder should be named after the provider, e.g. `panda`.