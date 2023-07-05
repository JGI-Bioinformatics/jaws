## Local Development

The new local development model uses docker images for all services, these images can be built/stored locally for local development, and they are also built/stored in the gitlab registry. Because the gitlab container registry is a new dependency, acquiring and configuring local access to the gitlab repository is needed.

* Gitlab Container Registry general info: https://docs.gitlab.com/ee/user/packages/
* Our Gitlab Registry: https://library.jgi.doe.gov/groups/advanced-analysis/-/packages
* Background info on docker registries: https://docs.docker.com/registry/introduction/

To build any images, you will need to have a Gitlab personal access token (or comparable) and use that to
authenticate your docker access.

* Getting a personal access token: https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html#create-a-personal-access-token
* Using the token to authenticate to gitlab's docker registry: https://docs.gitlab.com/ee/user/packages/container_registry/#authenticate-with-the-container-registry

This diagram (from [Microsoft](https://docs.microsoft.com/en-us/dotnet/architecture/microservices/docker-application-development-process/docker-app-development-workflow) ) shows the overall workflow for local development using Docker:

![Docker workflow](https://docs.microsoft.com/en-us/dotnet/architecture/microservices/docker-application-development-process/media/docker-app-development-workflow/life-cycle-containerized-apps-docker-cli.png)

You will need to do the following steps to use this workflow:

1. Make sure you have a functioning docker environment on your computer
2. Authenticate to the container registry on library.jgi.doe.gov
3. Use the build.sh script to build local copies of your images (this can result in pulling in dependencies from library.jgi.doe.gov)
   1. The site, central and dashboard images all depend on the rpc docker image, so it needs to be built before attempting to build site, central or dashboard
4. Use the docker-compose.yml file to start up a local environment using the docker images built in the previous step
5. Test your services locally
6. Make any changes required and loop back to step 3 to iterate until you have something you want to commit and push back to the repo
7. Push your changes back to the main repo - note that you should avoid pushing changes to docker-compose.yml which are specific to your setup, or which contains credentials

To deploy the JAWS services (except for the backend) locally you will want to make use
of the docker-compose.yml file and modify as needed. Some volume mounts are required in order for file uploads to
be processed by the containers. You will want to set your `DATA_HOME` with the locations of your
configuration files, log files and data upload files. Here is an example: 


```console
❯ pwd
/Users/mamelara/data/jaws
~/data/jaws on ☁️  (us-west-2)
❯ tree -L 2
.
├── configs
│   ├── auth
│   ├── central
│   ├── central-rpc
│   ├── cromwell
│   ├── daemon
│   └── site-central
├── laptop
│   ├── hello_world
│   └── output
├── local
│   ├── cromwell-workflow-logs
│   └── uploads
└── logs
    ├── auth
    ├── central
    ├── central-rpc
    ├── cromwell
    ├── cromwell-workflow-logs
    ├── daemon
    ├── jaws-auth
    └── site-central

22 directories, 0 files
> export DATA_HOME=/Users/mamelara/data/jaws
```

A sample configuration based on Mario's local test instance can be found in under /test/localdev of this repo. It can serve as a starting point for your own local test instance

The important thing to note is the name of the files, though these can be easily changed
in the `docker-compose.yml` file if you wish to use a different naming scheme. Once you've
created the following directory tree, the next step is to build the images.

The `laptop` and `local` directories are just there to "mock" locations on different filesystems. You can name these whatever directory
you want as long as you make sure the volume mounts match the names.

### Image building
At the root of the JAWS directory is a `build.sh` script that takes arguments for 
 the service name, the version and environment (eg - dev, prod). If you want to build central you will
want to run the following command: 

### ToDo: Explain how the $version param to build.sh must match version a version in the gitlab repo
### ToDo: Add details about authenticating access to the gitlab registry to download dependencies, or else make those images public

```console
> ./build.sh central 2.6.0 dev
```

This will build the development environment for your container script. This is useful if you want the code you are actively working on to be deployed in the container.  
Because rpc is a dependency for central, and site services, you will first want to build rpc docker container. 

```console
> ./build.sh rpc 2.6.0 dev
```

Once your images are built you can then modify the `docker-compose.yml` with the apppropriate 
image tags and then deploy the services. You'll want to deploy the mysql and rabbitmq services
first: 

```console
> docker-compose up -d db rabbitmq
```
You will want to setup your database first by creating some databases in mysql: 
You can enter the container either with `docker exec` or with the docker desktop. It is recommended 
to use the docker desktop since it provides a nice GUI.
```console
> docker exec -it <CONTAINER ID or CONTAINER NAMES> /bin/bash
```

If you're not sure of the container id or the container name you can run `docker ps`.
You will want to create the databases inside mysql: 

```console
$ mysql -u root -p 
$ Password: ********
$ CREATE DATABASE IF NOT EXISTS jaws_cromwell_dev;
$ CREATE DATABASE IF NOT EXISTS jaws_local_dev;
$ GRANT ALL PRIVILEGES ON * . * TO 'jaws'@'%';
```

Next you can bring up the other services without specifying their names.

```console
> docker-compose up -d
```

You should now have a deployed version of JAWS (except for the backend) on your local dev environment. To confirm
visit: http://localhost:80/api/v2/status and you should see the JAWS service statuses displayed in your browser. 

In order to use this local instance you will want to create a user in the mysql container. Enter the container and then
open up the mysql console. There is a table in the `jaws_central_dev` database called users. You will want to enter your
user information there. 

````console
$ mysql -u jaws -p
$ Password:
$ use jaws_central_dev;
$ insert into users (id, email, name, is_admin, is_dashboard, jaws_token) values ("yourusername", "youremail", true, false, "yourtoken");1
````

This should insert the first user in the database which will be your user. 

### Local Development Flow
The real power of deploying locally comes in being able to have your own instance of JAWS available. You have
all other services available for you with little setup. Let's look at an example on how to develop JAWS central code 
locally.

You'll want to use the `./build.sh` script to make sure you build the dev container version.

```console
./build.sh central 2.6.0 dev
```

This will build the central dev image with the tag `2.6.0-dev`. Next you will want to make sure you bind mount the source code you are working on: 

```yaml

volumes:
	- $DATA_HOME/configs/central/jaws-central.conf:/etc/config/jaws-central.conf
	- $DATA_HOME/laptop:/data/jaws/laptop:rw
	- $DATA_HOME/local:/data/jaws/local:rw
	- $DATA_HOME/logs/central:/var/log/jaws-central:rw
	- ./central:/usr/app
```

Note that the working directory is located in `/usr/app`. This is where the source code will be installed in the container.

You'll also want to just run an interactive shell without any commands running:

```yaml
central:
	image: central:2.6.0-dev
	stdin_open: true  # add these to the docker-compose file
	tty: true
```

Then you can run a `docker-compose up -d`. And your container will be running. You can enter the shell once again
with `docker exec -it <CONTAINER_ID> /bin/bash`. Once you are in there you can install your package inside
the container by running `pip install -e .`

Since you mounted your code, any changes reflected in your code will be reflected in the container. You can then
run that service with all the updated code changes you made to do some live testing. 
