# Jaws Site 

Jaws site is an RPC server that communicates via REST endpoints using AMQP
. It also connects with a SQL database to query and update for the status of
 each cromwell task. 
 
### Steps to install client
see the section "Installing all services" in the main README.md

## Starting the RPC server
### Configuration
Jaws site comes with a cli.py file that is the CLI client for the site
 . Before starting the server it is recommended you create a configuration
  file with the following contents: 
  
   - **AMQP**: Configuration for the RabbitMQ server. This includes the
    following parameters: 
        - host: the hostname where the rabbitmq server is located.
        - vhost: virtual host name.
        - user: username for the rabbitmq server.
        - password: password for the rabbitmq server.
        - queue: name of the messaging queue to send messages through.
   - **RPC**: configuration settings for the RPC server. 
        - host: hostname where site will be located.
        - vhost: virtual hostname.
        - password: password for the rpc server
        - queue: name of the rabbitmq queue
        - num_threads: number of consumers that will consume messages from
         queue.
        - max_retries: the maximum number of retries to connect rabbitmq
   - **GLOBUS**: Globus configuration settings. For more information check
    the Globus documentation.
        - client_id: ID for the client.
        - endpoint_id: name for the Globus endpoint. 
        - root_dir: root directory where data files where be placed.
        - subdir: name of the subdirectory where data files will be placed.
   - **DB**: Database configuration for SQL db.
        - host: name of the host where the database is located.
        - port: port where the database is listening to.
        - user: user for the database.
        - password: password for the database.
        - db: the name of the database.
        - dialect: dialect driver to connect to the database. For more
         information check out the sqlalchemy documentation.
   - **CROMWELL**: Configuration settings for Cromwell. 
        - host: name of the host where the cromwell service is running.
        - port: port to connect with cromwell
        - workflows_url: the workflows endpoint for using Cromwell's
         Workflows API. 
        - engine_status_url: endpoint where you can use Cromwell's Status API.
        

You can specify the location of a configuration file by setting the
 `JAWS_SITE_CONFIG` environment variable. 
 
### Setting up Jaws site services (RabbitMQ, Cromwell, MySQL)
Jaws comes with a `docker-compose.yml` that can be used to spin up services
via localhost. To use this file make sure you have `Docker` installed as
well as `docker-compose`. Once you have `docker-compose` installed, you can
spin up the services by being in the same directory as your `docker
-compose.yml` file and running the command `docker-compose up -d`. This
will start up the services as a background process. To confirm the status
of these services use `docker ps`.  

### Using the cli.py client 
Jaws comes with a `cli.py` file that can be used to manually test the Jaws
 site service. You will likely want to set PYTHONPATH to point to the root of
site.  

```shell script
export PYTHONPATH=/location/to/site/
```

Then you can start using the `cli.py` script.  

```shell script
Usage: cli.py [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  cancel-run     Cancel a run.
  check-runs     Check if any runs need action.
  daemon         Start JAWS-Site daemon.
  failure-logs   Get the logs for failed tasks.
  run-logs       Get the logs for a run.
  run-metadata   Get Cromwell metadata for a run.
  server         Start JAWS-Site RPC-server.
  server-status  Check Cromwell server status.
  submit-run     Submit a new run.
  task-ids       Retrieve the task IDs for a run.
  task-status    Get status of a run.
```

To start up the RPC server, run the command `./cli.py server`.  

