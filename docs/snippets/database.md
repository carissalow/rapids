??? info "Setting up a DATABASE_GROUP and its connection credentials."

    1. If you haven't done so, create an empty file called `#!bash credentials.yaml` in your RAPIDS root directory: 

    2. Add the following lines to `credentials.yaml` and replace your database-specific credentials (user, password, host, and database):

        ``` yaml
        MY_GROUP:
          database: MY_DATABASE
          host: MY_HOST
          password: MY_PASSWORD
          port: 3306
          user: MY_USER
        ```

    1. Notes
    
        1. The label `[MY_GROUP]` is arbitrary but it has to match the `[DATABASE_GROUP]` attribute of the data stream you choose to use.

        2. Indentation matters

        3. You can have more than one credentials group in `credentials.yaml`

    ??? hint "Upgrading from `./.env` from RAPIDS 0.x"
        In RAPIDS versions 0.x, database credentials were stored in a `./.env` file. If you are migrating from that type of file, you have two options:

        1. Migrate your credentials by hand:

            === "change .env format"
    
                ``` yaml
                [MY_GROUP]
                user=MY_USER
                password=MY_PASSWORD
                host=MY_HOST
                port=3306
                database=MY_DATABASE
                ```   
            
            === "to credentials.yaml format"
            
                ``` yaml
                MY_GROUP:
                  user: MY_USER
                  password: MY_PASSWORD
                  host: MY_HOST
                  port: 3306
                  database: MY_DATABASE
                ```

        2. Use the migration script we provide (make sure your conda environment is active):

            ```python
            python tools/update_format_env.py
            ```

    ??? hint "Connecting to localhost (host machine) from inside our docker container."
        If you are using RAPIDS' docker container and Docker-for-mac or Docker-for-Windows 18.03+, you can connect to a MySQL database in your host machine using `host.docker.internal` instead of `127.0.0.1` or `localhost`. In a Linux host, you need to run our docker container using `docker run --network="host" -d moshiresearch/rapids:latest` and then `127.0.0.1` will point to your host machine.