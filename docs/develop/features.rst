Add new features to RAPIDS
============================

Take accelerometer features as an example.

#. Add your script to accelerometer_ folder

    - Copy the signature of the base_accelerometer_features() function_ for your own feature function

#. Add any parameters you need for your function

    - Add your parameters to the settings_ of accelerometer sensor in config file
    - Add your parameters to the params_ of accelerometer_features rule in features.snakefile

#. Merge your new features with the existent features

    - Call the function you just created below this line (LINK_) of accelerometer_features.py script

#. Update config file

    - Add your new feature names to the ``FEATURES`` list for accelerometer in the config_ file

.. _accelerometer: https://github.com/carissalow/rapids/tree/master/src/features/accelerometer
.. _function: https://github.com/carissalow/rapids/blob/master/src/features/accelerometer/accelerometer_base.py#L35
.. _settings: https://github.com/carissalow/rapids/blob/master/config.yaml#L100
.. _params: https://github.com/carissalow/rapids/blob/master/rules/features.snakefile#L146
.. _LINK: https://github.com/carissalow/rapids/blob/master/src/features/accelerometer_features.py#L10
.. _config: https://github.com/carissalow/rapids/blob/master/config.yaml#L102
