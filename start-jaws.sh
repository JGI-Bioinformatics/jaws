#!/bin/bash

docker-compose up -d central_db central_rabbitmq site_db site_rabbitmq cromwell

sleep 10  # sleep and wait for dbs and rabbitmqs to start-up.

docker-compose up -d
