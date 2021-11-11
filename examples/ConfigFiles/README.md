# Call-caching Explained
Cromwell can run either in server mode or non-server (client) mode. 

To use the call-cache option, you need to run cromwell with a database to keep track of the previous cromwell runs and there are two types of databases:
1. in-memory: this db will only be persistent as long as the cromwell server is running. If you run in client mode, there will be no call-caching.
2. mysql or sqlite: A persistent database will allow you to run cromwell either in server or client mode and still have call-caching. 

## Server Mode
An in-memory database comes with cromwell by default, which means if you run in server mode, you don't need to set up the database. 
You still need a cromwell.conf file to change the default settings to use call-cache.

### Submitting to a cromwell server 
```
# start the server with call-caching turned on
java -Dconfig.file=cromwell.conf -jar <path_to_your>/cromwell.jar server &

# submit the WDL to Cromwell. Call-caching is turned on via config file.
curl -X POST --header "Accept: application/json" -v "localhost:50011/api/workflows/v1" -F workflowSource=@hello.wdl -F workflowInputs=@inputs.json -F workflowOptions=@options.json

# submit with call-caching turned off
curl -X POST --header "Accept: application/json" -v "localhost:50011/api/workflows/v1" -F workflowSource=@hello.wdl --F workflowOptions=@options.json

# turn off call-caching by submitting with
-F workflowOptions=@options.json

# see the metadata to check if call caching was on or off
# use fq to parse the json
curl --header "Accept: application/json" -v "localhost:8088/api/workflows/v1/${wid}/metadata" | jq '.calls."test.hello"[0].callCaching'
```

To get a pretty json output, you can also pipe the json to 
```
python -m json.tool
or
json_pp -json_opt pretty,canonical

# for example
curl --header "Accept: application/json" -v "localhost:8088/api/workflows/v1/${wid}/metadata" | python -m json.tool
```

## Client Mode
For client mode, you need to set up a persistant database:  mysql, or in this example, a file-based db. You can configure these databases in a cromwell config which is then part of the submit command. The config examples supplied here (except parsl.conf) use the file-based 'HsqldbProfile' db.  See the cromwell docs for [HsqldbProfile](https://cromwell.readthedocs.io/en/stable/Configuring/) under the section "Using Cromwell with file-based database (No server required)".

```
java -Dconfig.file=cromwell_docker.conf \
-Dbackend.providers.local.config.dockerRoot=$(pwd)/cromwell-execution \
-Dbackend.providers.local.config.root=$(pwd)/cromwell-execution \
-jar cromwell.jar \
run hello.wdl -i hello.json
```

Note that you need to overwrite `dockerRoot` and `root` from the config file, to point to your own paths.

## The options.json File
You can overwrite certain config settings with the options.json file.  For example, you can turn off call-caching or use another backend.
For example, you can set "backend": [local|slurm].

You can also set the options on the command line with `-D` like 
```
java -Dbackend=slurm ...etc.
```
