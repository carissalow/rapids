# Configuration Schema
The configuration schema is a JSON schema, which is for validating the strcuture of sensor data. By running the command `snakemake --list-params-changes`, rapids will check if all settings match the standard without computing the sensor data. The configuration can be divided into three parts: required, definitions and properties. 

### Required
This part specifies the required features that need to be defined under the properties. If there is any missing feature, it will throw an error message. 
For example, there is a required feature called `Testing_feature` but it is not defined in the properties. When running `snakemake --list-params-changes`, here is the error message you will see, `ValidationError: 'Testing_feature' is a required property`. 

### Definitions 
Features might share some common properties, to reuse the structure of the common properties in other places, you can define the property under the definitions, then using `$ref` to refer to the defined property. 

For example, there is a defined property called `PROVIDERS` in `definitions` and a new feature called `phone_messages` has such property. 

![picture](/img/provider_example.png)

To reuse the definiation of `PROVIDERS`, all you need to do is adding `$ref: "#/definitions/PROVIDER"` to the additional property of `phone_messages`.

![picture](/img/phone_messages_example.png)

You can also overwrite the existing property in `PROVIDERS` if the definiation of a certain property is different from the current `PROVIDER`, or extend the properties of `PROVIDERS`.


For more examples and details, please refer to http://json-schema.org/understanding-json-schema/structuring.html#reuse

### Properties
This is where features are defined. The definition of a feature usually contains three things. The type of the feature, a list of required properties, and the definition of each property. A feature is usually an object, which is like the dictionary in Python. Everthing inside the object is a pair, meaning that it has a key and a value. In JSON schema, the `keys` must always be string. Besides object, other common types are string, numeric, array, boolean and null. 

For example,

![picture](/img/feature_example.png)

The sensor feature PHONE_AWARE_LOG is an object with two required properties, `TABLE` and `PROVIDERS`. The `TABLE` contains the name of the data table used for computation, which is a string. The `PROVIDERS` has two types, `null` and `object` meaning that PHONE_AWARE_LOG may or may not have a provider. 

### Refrences

- Understanding JSON Schema: http://json-schema.org/understanding-json-schema/index.html
- A Media Type for Describing JSON Documents: https://tools.ietf.org/html/draft-handrews-json-schema-01