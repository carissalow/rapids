# Validation schema of `config.yaml`

!!! hint "Why do we need to validate the `config.yaml`?"
    Most of the key/values in the `config.yaml` are constrained to a set of possible values or types. For example `[TIME_SEGMENTS][TYPE]` can only be one of `["FREQUENCY", "PERIODIC", "EVENT"]`, and `[TIMEZONE]` has to be a string. 
    
    We should show the user an error if that's not the case. We could validate this in Python or R but since we reuse scripts and keys in multiple places, tracking these validations can be time consuming and get out of control. Thus, we do these validations through a schema and check that schema before RAPIDS starts processing any data so the user can see the error right away.

    Keep in mind these validations can only cover certain base cases. Some validations that require more complex logic should still be done in the respective script. For example, we can check that a CSV file path actually ends in `.csv` but we can only check that the file actually exists in a Python script.
 
The structure and values of the `config.yaml` file are validated using a YAML schema stored in `tools/config.schema.yaml`. Each key in `config.yaml`, for example `PIDS`, has a corresponding entry in the schema where we can validate its type, possible values, required properties, min and max values, among other things. 

The `config.yaml` is validated against the schema every time RAPIDS runs (see the top of the `Snakefile`):

```python
validate(config, "tools/config.schema.yaml")
```

## Structure of the schema

The schema has three main sections `required`, `definitions`, and `properties`. All of them are just nested key/value YAML pairs, where the value can be a primitive type (`integer`, `string`, `boolean`, `number`) or can be another key/value pair (`object`).

### required
`required` lists `properties` that should be present in the `config.yaml`. We will almost always add every `config.yaml` key to this list (meaning that the user cannot delete any of those keys like `TIMEZONE` or `PIDS`). 

### definitions
`definitions` lists key/values that are common to different `properties` so we can reuse them. You can define a key/value under `definitions` and use `$ref` to refer to it in any `property`. 

For example, every sensor like `[PHONE_ACCELEROMETER]` has one or more providers like `RAPIDS` and `PANDA`, these providers have some common properties like the `COMPUTE` flag or the `SRC_SCRIPT` string. Therefore we define a shared provider "template" that is used by every provider and extended with properties exclusive to each one of them. For example:

=== "provider definition (template)"
    The `PROVIDER` definition will be used later on different `properties`.

    ```yaml
    PROVIDER:
        type: object
        required: [COMPUTE, SRC_SCRIPT, FEATURES]
        properties:
        COMPUTE:
            type: boolean
        FEATURES:
            type: [array, object]
        SRC_SCRIPT:
            type: string
            pattern: "^.*\\.(py|R)$"
    ```

=== "provider reusing and extending the template"
    Notice that `RAPIDS` (a provider) uses and extends the `PROVIDER` template in this example. The `FEATURES` key is overriding the `FEATURES` key from the `#/definitions/PROVIDER` template but is keeping the validation for `COMPUTE`, and `SRC_SCRIPT`. For more details about reusing properties, go to this [link](http://json-schema.org/understanding-json-schema/structuring.html#reuse)

    ```yaml hl_lines="9 10"
    PHONE_ACCELEROMETER:
        type: object
         # .. other properties
        PROVIDERS:
            type: ["null", object]
            properties:
            RAPIDS:
                allOf:
                - $ref: "#/definitions/PROVIDER"
                - properties:
                    FEATURES: 
                        type: array
                        uniqueItems: True
                        items:
                        type: string
                        enum: ["maxmagnitude", "minmagnitude", "avgmagnitude", "medianmagnitude", "stdmagnitude"]
    ```



### properties

`properties` are nested key/values that describe the different components of our `config.yaml` file. Values can be of one or more primitive types like `string`, `number`, `array`, `boolean` and `null`. Values can also be another key/value pair (of type `object`) that are similar to a dictionary in Python.

For example, the following property validates the `PIDS` of our `config.yaml`. It checks that `PIDS` is an `array` with unique items of type `string`.

```yaml
PIDS:
    type: array
    uniqueItems: True
    items:
      type: string
```

## Modifying the schema

!!! hint "Validating the `config.yaml` during development"
    If you updated the schema and want to check the `config.yaml` is compliant, you can run the command `snakemake --list-params-changes`. You will see `Building DAG of jobs...` if there are no problems or an error message otherwise (try setting any `COMPUTE` flag to a string like `test` instead of `False/True`).
    
    You can use this command without having to configure RAPIDS to process any participants or sensors.

You can validate different aspects of each key/value in our `config.yaml` file:

=== "number/integer"
    Including min and max values
    ```yaml
    MINUTE_RATIO_THRESHOLD_FOR_VALID_YIELDED_HOURS:
        type: number
        minimum: 0
        maximum: 1

    FUSED_RESAMPLED_CONSECUTIVE_THRESHOLD:
        type: integer
        exclusiveMinimum: 0
    ```
=== "string"
    Including valid values (`enum`)
    ```yaml
    items:
        type: string
        enum: ["count", "maxlux", "minlux", "avglux", "medianlux", "stdlux"]
    ```
=== "boolean"
    ```yaml
    MINUTES_DATA_USED:
        type: boolean
    ```
=== "array"
    Including whether or not it should have unique values, the type of the array's elements (`strings`, `numbers`) and valid values (`enum`).
    ```yaml
    MESSAGES_TYPES:
        type: array
        uniqueItems: True
        items:
            type: string
            enum: ["received", "sent"]
    ```
=== "object"
    `PARENT` is an object that has two properties. `KID1` is one of those properties that are, in turn, another object that will reuse the  `"#/definitions/PROVIDER"` `definition` **AND** also include (extend) two extra properties `GRAND_KID1` of type `array` and `GRAND_KID2` of type `number`. `KID2` is another property of `PARENT` of type `boolean`.

    The schema validation looks like this
    ```yaml
    PARENT:
        type: object
        properties:
          KID1:
            allOf:
              - $ref: "#/definitions/PROVIDER"
              - properties:
                  GRAND_KID1:
                    type: array
                    uniqueItems: True
                  GRAND_KID2:
                    type: number
          KID2:
            type: boolean
    ```

    The `config.yaml` key that the previous schema validates looks like this:
    ```yaml
    PARENT:
        KID1:
            # These four come from the `PROVIDER` definition (template)
            COMPUTE: False
            FEATURES: [x, y] # an array
            SRC_SCRIPT: "a path to a py or R script"

            # This two come from the extension
            GRAND_KID1: [a, b] # an array
            GRAND_KID2: 5.1 # an number
         KID2: True # a boolean
    ```

## Verifying the schema is correct
We recommend that before you start modifying the schema you modify the `config.yaml` key that you want to validate with an invalid value. For example, if you want to validate that `COMPUTE` is boolean, you set `COMPUTE: 123`. Then create your validation, run `snakemake --list-params-changes` and make sure your validation fails (123 is not `boolean`), and then set the key to the correct value. In other words, make sure it's broken first so that you know that your validation works.

!!! warning
    **Be careful**. You can check that the schema `config.schema.yaml` has a valid format by running `python tools/check_schema.py`. You will see this message if its structure is correct: `Schema is OK`. However, we don't have a way to detect typos, for example `allOf` will work but `allOF` won't (capital `F`) and it won't show any error. That's why we recommend to start with an invalid key/value in your `config.yaml` so that you can be sure the schema validation finds the problem.

## Useful resources

Read the following links to learn more about what we can validate with schemas. They are based on `JSON` instead of `YAML` schemas but the same concepts apply.

- [Understanding JSON Schemas](http://json-schema.org/understanding-json-schema/index.html)
- [Specification of the JSON schema we use](https://tools.ietf.org/html/draft-handrews-json-schema-01)